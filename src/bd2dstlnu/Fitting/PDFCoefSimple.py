import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import json

class PDFAngularCoef:
    def __init__(self):
        self._data = None
        self._Ndata = 0.0
        self._angint = None
        self._dataErr = None
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
    def angint(self) -> np.ndarray:
        """Angular integrals correpsonding to the data's binning scheme"""
        return self._angint
    @angint.setter
    def angint(self, a: np.ndarray):
        self._angint = a
    @property
    def m(self) -> Minuit:
        """Minuit minimizer object"""
        return self._m

    def setData(self, data: np.ndarray, data_err: np.ndarray = None):
        self._data = data
        self._Ndata = np.sum(data)
        if type(data_err) != type(None):
            self._dataErr = data_err
        else:
            self._dataErr = np.where(data > 0, np.sqrt(data), 1.0)

    def pdf(self, params: list) -> float:
        # params in the order ['1s', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]
        # Set 1c by normalization
        I1c = (4 - 6*params[0] + params[2] + 2*params[1])/3
        j_vec = np.concatenate(([params[0], I1c], params[1:]))
        hist = np.dot(self.angint, j_vec)*(9/(32*np.pi))
        return hist
    
    def chi2(self, params: list) -> float:
        # Simple chi2 fit - assume uncorrletad bin yields but this will likely not be the case
        pred = self.pdf(params)
        h_pred = self.Ndata*(pred/np.sum(pred))
        chiSq = np.sum((self.data - h_pred)**2/(self.dataErr**2))
        return chiSq
    
    def fit(self, params_initial: list[float], params_limits: list[float] = []):

        m = Minuit(self.chi2, params_initial)
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
    
    def saveResults(self, savepath: str, params: list[float], cov: np.ndarray, optionals: dict = {}):
        res = {"val" : params, "cov" : cov.tolist(), "fmin" : self.m.fval, "valid" : self.m.valid, "covQual" : self.m.accurate, **optionals}
        with open(savepath, "w") as outfile:
            json.dump(res, outfile, indent=2)
        return res