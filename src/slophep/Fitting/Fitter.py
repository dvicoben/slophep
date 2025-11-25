import numpy as np
from iminuit import Minuit
from slophep.Fitting.Parameter import Parameter, ParameterManager
from slophep.Fitting.CostFuncs import FigureOfMerit
from seaborn import heatmap
import matplotlib.pyplot as plt

class imFitter:
    def __init__(self, 
                 FoM: FigureOfMerit):
        
        self._fom = FoM
        self._m = None

    @property
    def fom(self) -> FigureOfMerit: 
        """Figure of merit to minimise"""
        return self._fom
    @property
    def param_manager(self) -> ParameterManager: 
        """Parameter manager"""
        return self.fom.param_manager
    @property
    def m(self) -> Minuit: 
        """iminuit object"""
        return self._m
    @property
    def param_list(self) -> list[str]: 
        """Parameter list"""
        return self.fom.param_list

    def migrad(self):
        """Minimise FoM using migrad"""
        vals = [self.param_manager.getVal(ipar) for ipar in self.param_list]
        self._m = Minuit(self.fom, vals)
        for i in range(len(self.param_list)):
            ipar = self.param_manager.getParam(self.param_list[i])
            self.m.limits[i] = (ipar.limlo, ipar.limhi)
            if ipar.fixed:
                self.m.fixto(i, ipar.val)
        
        ncall = 1000000
        self.m.migrad(ncall=ncall)

    def hesse(self):
        """Calculate Hessian"""
        self.m.hesse()

    def fit(self, sumW2error: bool = False):
        """Perform fit minimisation and calculate hessian"""
        self.migrad()
        self.hesse()
        
        res_list = [elem.value for elem in self.m.params]
        err = [elem.error for elem in self.m.params]
        if sumW2error:
            m2 = Minuit(self.fom.calc_nllW2, res_list)
            m2.hesse()
            c_inv = np.linalg.inv(m2.covariance)
            cov = np.matmul(self.m.covariance, np.matmul(c_inv, self.m.covariance))
            err = np.sqrt(np.diag(cov))


        print("Function Minimum:", self.m.fval)
        if self.m.valid:
            print("Fit converged")
        else:
            print("Fit did not converge")

        if self.m.accurate:
            print("Covariance matrix accurate")
        else:
            print("Covariance matrix not accurate")
        
        return res_list, err
    
    def corr(self, savepath: str = None, labels: list[str] = []) -> np.ndarray:
        """Parameter correlation matrix

        Parameters
        ----------
        savepath : str, optional
            Path to save correlation matrix to, if desired, by default None
        labels : list[str], optional
            Labels to use for parameters if not their names, by default []

        Returns
        -------
        np.ndarray
            Correlation matrix
        """
        corr = self.m.covariance.correlation()
        if savepath:
            fig, ax = plt.subplots(1, 1)
            labs = labels if len(labels) > 0 else self.param_list
            heatmap(data=corr, vmin=-1, vmax=1, annot=True, cmap='bwr', fmt=".3f", ax=ax,
                    xticklabels=labs, yticklabels=labs)
            ax.tick_params(axis='x', labelrotation=90)
            fig.savefig(savepath, bbox_inches="tight", dpi=100)
            plt.close()
        return corr