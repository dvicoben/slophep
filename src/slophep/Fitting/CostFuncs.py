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
    def param_manager(self) -> ParameterManager: 
        """The parameter manager"""
        return self._param_manager
    @property
    def param_list(self) -> list[str]: 
        """List of parameters"""
        return self._param_list
    @property
    def pdf(self) -> PDFBase: 
        """PDF to calculate FoM for"""
        return self._pdf
    @property
    def data(self) -> np.ndarray: 
        """Data to fit"""
        return self._data
    @property
    def dataNorm(self) -> float: 
        """Sum of data yield"""
        return self._dataNorm
    @property
    def dataErr(self) -> np.ndarray: 
        """Uncertainties of data"""
        return self._dataErr
    @property
    def scaleFactorW2(self) -> np.ndarray: 
        """Scale factor for sumW2 likelihood"""
        return self._scaleFactorW2
    @property
    def ignore_constraint(self) -> bool: 
        """Whether to ignore constraint terms in the likelihood/figure of merit"""
        return self._ignore_constraint

    def setVal(self, parname_man: str, val: float):
        """Set parameter value"""
        self.param_manager.setVal(parname_man, val)

    def setVals(self, pars_man: dict[str, float]):
        """Set value for several parameters, with dictionary {manager_name : value}"""
        self.param_manager.setVals(pars_man)

    def updateVals(self, paramvals: list[float]):
        """Update all parameter values with list input in order of self.param_list"""
        params = {parname : parval for parname, parval in zip(self.param_list, paramvals)}
        self.param_manager.setVals(params)

    def get_model(self) -> np.ndarray:
        """Fit model"""
        if self._normPDF_before_eval:
            model = self.pdf.normpdf()
            model_scaled = self.dataNorm*model
            return model_scaled
        return self.pdf.pdf()

    def get_constraint_term(self) -> float:
        """Get FoM constraint term

        Returns
        -------
        float
            Gaussian likelihood term. Return 0 if ignore\_constraint == True.
        """
        if self.ignore_constraint:
            return 0.0
        
        consterm = 0.0
        for iparcode, ipar in self.pdf.param_manager.params.items():
            if ipar.constrain:
                consterm += cf.nll_gaussian(ipar.getVal(), ipar.constraintvals[0], ipar.constraintvals[1])
        return consterm

    def __call__(self, paramvals: list[float]) -> float:
        """Update parameters and evaluate FoM

        Parameters
        ----------
        paramvals : list[float]
            Parameter values to evaluate at, in order of param_list

        Returns
        -------
        float
            Figure of merit

        Raises
        ------
        Exception
            Figure of merit call implemented in child class
        """
        self.updateVals(paramvals)
        # return self.nll()
        raise Exception("FOM __call__ specified in derived class")

    def nll(self) -> float:
        """Returns -2LL"""
        return cf.poisson_chi2(self.data, self.get_model()) + self.get_constraint_term()
    
    def nllW2(self) -> float:
        """Return -2LL with SumW2, to use for uncertainties"""
        return cf.poisson_chi2_W2(self.data, self.get_model(), self.scaleFactorW2) + self.get_constraint_term()
    
    def calc_nllW2(self, paramvals: list[float]) -> float:
        """Calculate -2LL with SumW2 for particular set of parameters

        Parameters
        ----------
        paramvals : list[float]

        Returns
        -------
        float
        """
        self.updateVals(paramvals)
        return self.nllW2()

    def chi2_uncorrelated(self) -> float:
        """Uncorrelated chi2"""
        chi2 = np.sum((self.data - self.get_model())**2 / (self.dataErr**2))
        return chi2 + self.get_constraint_term()
    


class Chi2Uncorrelated(FigureOfMerit):
    def __init__(self,
                 pdf: PDFBase,
                 param_manager: ParameterManager,
                 data: np.ndarray,
                 dataErr: np.ndarray,
                 normPDF_before_eval: bool = False,
                 ignore_constraint: bool = False):
        if np.shape(dataErr) != np.shape(data):
            raise ValueError(f"Data and errors do not have the same shape {np.shape(data)} {np.shape(dataErr)}")
        super().__init__(pdf=pdf, param_manager=param_manager, data=data, dataErr=dataErr, 
                         normPDF_before_eval=normPDF_before_eval, ignore_constraint=ignore_constraint)

    def __call__(self, paramvals: list[float]) -> float:
        """The uncorrelated chi2"""
        self.updateVals(paramvals)
        return self.chi2_uncorrelated()



class Chi2Correlated(FigureOfMerit):
    def __init__(self,
                 pdf: PDFBase,
                 param_manager: ParameterManager,
                 data: np.ndarray,
                 cov: np.ndarray,
                 normPDF_before_eval: bool = False,
                 ignore_constraint: bool = False):
        if len(np.shape(data)) > 1:
            raise ValueError(f"Data has shape {np.shape(data)} - must be flattened to 1D")
        if np.shape(cov) != (np.shape(data), np.shape(data)):
            raise ValueError(f"Cov matrix shape {np.shape(cov)}, expected {(np.shape(data), np.shape(data))}")
        
        super().__init__(pdf=pdf, param_manager=param_manager, data=data, dataErr=np.sqrt(np.diag(cov)), 
                         normPDF_before_eval=normPDF_before_eval, ignore_constraint=ignore_constraint)
        self._cov = cov
        self._invcov = np.linalg.inv(self.cov)
    
    @property
    def cov(self) -> np.ndarray:
        """Covariance matrix"""
        return self._cov
    @property
    def invcov(self) -> np.ndarray:
        """Inverse covariance matrix"""
        return self._invcov

    def chi2_correlated(self) -> float:
        """Return chi2 for correlated bins, (obs - pred) * invcov * (obs - pred)"""
        delta_p = (self.data - self.get_model())
        chi2 = np.dot(delta_p, np.dot(self.invcov, delta_p))
        return chi2 + self.get_constraint_term()

    def __call__(self, paramvals: list[float]) -> float:
        """Correlated chi2"""
        self.updateVals(paramvals)
        return self.chi2_correlated()