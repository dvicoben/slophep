import numpy as np

import slophep.Fitting.CostFuncMath as cf
from slophep.Fitting.Parameter import Parameter, ParameterManager
from slophep.Fitting.PDFBase import PDFBase

class FigureOfMerit:
    def __init__(self, 
                 pdf: PDFBase,
                 param_manager: ParameterManager,
                 data: np.ndarray,
                 dataErr: np.ndarray = np.array([]),
                 normPDF_before_eval: bool = False,
                 ignore_constraint: bool = False):
        
        self._data = data
        self._dataNorm = np.sum(self.data)
        self._dataErr = dataErr
        self._scaleFactorW2 = (dataErr**2)/np.where(data != 0, data, 1.0)
        self._ignore_constraint = ignore_constraint

        self._param_manager = param_manager
        self._param_list = [ikey for ikey in self.param_manager.id_map]
    
        self._pdf = pdf
        self._pdf.set_param_manager(self._param_manager)
        self._pdf.link_params_check()
        self._normPDF_before_eval = normPDF_before_eval

    @property
    def param_manager(self) -> ParameterManager: return self._param_manager
    @property
    def param_list(self) -> list[str]: return self._param_list
    @property
    def pdf(self) -> PDFBase: return self._pdf
    @property
    def data(self) -> np.ndarray: return self._data
    @property
    def dataNorm(self) -> float: return self._dataNorm
    @property
    def dataErr(self) -> np.ndarray: return self._dataErr
    @property
    def scaleFactorW2(self) -> np.ndarray: return self._scaleFactorW2
    @property
    def ignore_constraint(self) -> bool: return self._ignore_constraint

    def setVal(self, parname_man: str, val: float):
        self.param_manager.setVal(parname_man, val)

    def setVals(self, pars_man: dict[str, float]):
        self.param_manager.setVals(pars_man)

    def updateVals(self, paramvals: list[float]):
        params = {parname : parval for parname, parval in zip(self.param_list, paramvals)}
        self.param_manager.setVals(params)

    def get_model(self):
        if self._normPDF_before_eval:
            model = self.pdf.normpdf()
            model_scaled = self.dataNorm*model
            return model_scaled
        return self.pdf.pdf()

    def get_constraint_term(self):
        if self.ignore_constraint:
            return 0.0
        
        consterm = 0.0
        for iparcode, ipar in self.pdf.param_manager.params.items():
            if ipar.constrain:
                consterm += cf.nll_gaussian(ipar.getVal(), ipar.constraintvals[0], ipar.constraintvals[1])
        return consterm

    def __call__(self, paramvals: list[float]) -> float:
        self.updateVals(paramvals)
        # return self.nll()
        raise Exception("FOM __call__ specified in derived class")

    def nll(self):
        return cf.poisson_chi2(self.data, self.get_model()) + self.get_constraint_term()
    
    def nllW2(self):
        return cf.poisson_chi2_W2(self.data, self.get_model(), self.scaleFactorW2) + self.get_constraint_term()
    
    def calc_nllW2(self, paramvals: list[float]):
        self.updateVals(paramvals)
        return self.nllW2()

    def chi2_uncorrelated(self):
        chi2 = np.sum((self.data - self.get_model())**2 / (self.dataErr**2))
        return chi2 + self.get_constraint_term()
    


class Chi2Uncorrelated(FigureOfMerit):
    def __init__(self,
                 pdf: PDFBase,
                 param_manager: ParameterManager,
                 data: np.ndarray,
                 dataErr: np.ndarray,
                 normPDF_before_eval: bool = False):
        if np.shape(dataErr) != np.shape(data):
            raise ValueError(f"Data and errors do not have the same shape {np.shape(data)} {np.shape(dataErr)}")
        super().__init__(pdf=pdf, param_manager=param_manager, data=data, dataErr=dataErr, normPDF_before_eval=normPDF_before_eval)

    def __call__(self, paramvals: list[float]) -> float:
        self.updateVals(paramvals)
        return self.chi2_uncorrelated()



class Chi2Correlated(FigureOfMerit):
    def __init__(self,
                 pdf: PDFBase,
                 param_manager: ParameterManager,
                 data: np.ndarray,
                 cov: np.ndarray,
                 normPDF_before_eval: bool = False):
        if len(np.shape(data)) > 1:
            raise ValueError(f"Data has shape {np.shape(data)} - must be flattened to 1D")
        if np.shape(cov) != (np.shape(data), np.shape(data)):
            raise ValueError(f"Cov matrix shape {np.shape(cov)}, expected {(np.shape(data), np.shape(data))}")
        
        super().__init__(pdf=pdf, param_manager=param_manager, data=data, dataErr=cov, normPDF_before_eval=normPDF_before_eval)
        self._invcov = np.linalg.inv(self.dataErr)
    
    @property
    def invcov(self) -> np.ndarray:
        return self._invcov

    def chi2_correlated(self):
        delta_p = (self.data - self.get_model())
        chi2 = np.dot(delta_p, np.dot(self.invcov, delta_p))
        return chi2 + self.get_constraint_term()

    def __call__(self, paramvals: list[float]) -> float:
        self.updateVals(paramvals)
        return self.chi2_correlated()