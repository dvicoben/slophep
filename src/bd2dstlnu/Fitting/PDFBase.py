import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import json
import pickle


class PDFBase:
    def __init__(self):
        self._data = None
        self._dataErr = None
        self._dataCov = None
        self._invCov = None
        self._Ndata = 0.0
        self._m = None

    @property
    def data(self) -> np.ndarray:
        """Data in an NDarray in ctd, ctl, chi"""
        return self._data
    @property
    def Ndata(self) -> np.ndarray:
        """Sum of all bins, to use for normalisation, cached to avoid recomputing"""
        return self._Ndata
    @property
    def dataErr(self) -> np.ndarray:
        """Data uncertainty"""
        return self._dataErr
    @property
    def dataCov(self) -> np.ndarray:
        """Data covariance"""
        return self._dataCov
    @property
    def invCov(self) -> np.ndarray:
        """Inverse of data covariance matrix"""
        return self._invCov
    @property
    def m(self) -> Minuit:
        """Minuit minimizer object"""
        return self._m
    
    def setData(self, data: np.ndarray, data_err: np.ndarray = None, errtype: str = "std"):
        self._data = data
        self._Ndata = np.sum(data)
        if errtype == "std":
            if type(data_err) != type(None):
                self._dataErr = data_err
            else:
                self._dataErr = np.where(data > 0, np.sqrt(data), 1.0)
        elif errtype == "cov":
            self._dataCov = data_err
            self._invCov = np.linalg.inv(data_err)

    def pdf(self, params: list) -> float:
        raise NotImplementedError("PDF should be implemented in derived class")

    def chi2(self, params: list) -> float:
        pred = self.pdf(params)
        h_pred = self.Ndata*(pred/np.sum(pred))
        delta_p = (self.data - h_pred) # may need to unroll/flatten histograms
        chiSq = np.dot(delta_p, np.dot(self._invCov, delta_p))
        return chiSq

    def chi2_uncorrelated(self, params: list) -> float:
        # Simple chi2 fit - assume uncorrletad bin yields but this will likely not be the case
        pred = self.pdf(params)
        h_pred = self.Ndata*(pred/np.sum(pred))
        chiSq = np.sum((self.data - h_pred)**2/(self.dataErr**2))
        return chiSq
    
    def fit(self, params_initial: list[float], params_limits: list[float] = [], fom: str = "chi2_uncorrelated"):

        m = Minuit(getattr(self, fom), params_initial)
        for i in range(len(params_limits)):
            if params_limits[i] == None: 
                continue
            else: 
                m.limits[i] = params_limits[i]
        
        m.migrad()
        m.hesse()
        print("Function Minimum:", m.fval)
        if m.valid:
            print("Fit converged")
        else:
            print("Fit did not converge")

        if m.accurate:
            print("Covariance matrix accurate")
        else:
            print("Covariance matrix not accurate")

        self._m = m
        res_list = [elem.value for elem in m.params]
        cov = np.array(m.covariance)
        return res_list, cov
    
    def saveResults(self, savepath: str, optionals: dict = {}):
        if type(self.m) == type(None):
            print("WARNING: No minimization, no results to save!")
            return {}

        params = [elem.value for elem in self.m.params]
        error = [elem.error for elem in self.m.params]
        cov = np.array(self.m.covariance)
        res = {
            "val" : params, 
            "err": error, 
            "cov" : cov.tolist(), 
            "fmin" : self.m.fval, 
            "valid" : self.m.valid, 
            "covQual" : self.m.accurate, 
            **optionals
        }
        if not savepath:
            return res

        if ".json" in savepath:
            with open(savepath, "w") as outfile:
                json.dump(res, outfile, indent=2)
        else:
            with open(savepath, "wb") as outfile:
                    pickle.dump(res, outfile)
        
        return res