import numpy as np
from slophep.Fitting.Parameter import Parameter, ParameterManager
from slophep.Predictions.Observables import ObservableBase

class PDFBase:
    def __init__(self, 
                 name: str,
                 predobj: ObservableBase):
        
        self._name = name
        self._paramslink = {}
        self._param_manager = None
        self._params_list = []
        self._binning = None
        self._predobj = predobj

    @property
    def name(self) -> str: return self._name
    @property
    def param_manager(self) -> ParameterManager: return self._param_manager
    @property
    def params_list(self) -> list[str]: return self._params_list
    @property
    def params_link(self) -> dict: return self._paramslink
    @property
    def binning(self) -> np.ndarray: return self._binning
    @property
    def predobj(self) -> ObservableBase: 
        return self._predobj

    def set_param_manager(self, param_manager: ParameterManager):
        self._param_manager = param_manager

    def link_param(self, name_int: str, name_manager: str):
        if name_int not in self.params_list:
            raise KeyError(f"{name_int} is not a PDF parameter")
        # if name_manager not in self._param_manager.id_map:
        #     raise KeyError(f"{name_manager} is not in parameter_manager")
        self._paramslink[name_int] = name_manager
    
    def link_params(self, link_dict: dict = {}):
        if len(link_dict) < 1:
            for ipar in self.params_list:
                self.link_param(ipar, ipar)
            return
        
        for ipar in link_dict:
            self.link_param(ipar, link_dict[ipar])

    def link_params_check(self):
        print("Checking parmeter links")
        for ipar_int, ipar_man in self.params_link.items():
            if ipar_man not in self.param_manager.id_map:
                raise KeyError(f"{ipar_man} (linked to {ipar_int}) is not in parameter_manager")
            print(f"- internal: {ipar_int}, external: {ipar_man}, manager_id: {self.param_manager.id_map[ipar_man]}")
        print("Check Done!")

    def set_binning(self, binning: np.ndarray):
        self._binning = binning

    def getVal(self, parname: str):
        parname_man = self.params_link[parname]
        return self.param_manager.getVal(parname_man)
    
    def setVal(self, parname: str, val: float):
        parname_man = self.params_link[parname]
        self.param_manager.setVal(parname_man, val)

    def setVals(self, pars: dict):
        pars_man = {self.params_link[k] : pars[k] for k in pars}
        self.param_manager.setVals(pars_man)

    def pdf(self) -> np.ndarray:
        raise NotImplementedError("PDF implemented in child class!")

    def normpdf(self) -> np.ndarray:
        pdf = self.pdf()
        return pdf/np.sum(pdf)
