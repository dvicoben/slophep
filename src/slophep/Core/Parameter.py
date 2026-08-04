from __future__ import annotations
import flavio
from typing import TypeVar, Any

from slophep.Core.physdata import DEFAULT_PARAMS
import copy
import logging

logger = logging.getLogger(__name__)

T = TypeVar("T")

class Parameter:
    def __init__(self, 
                 name: str,
                 val: T):
        """_summary_

        Parameters
        ----------
        name : str
            Parameter name
        val : T
            Parameter value, arbitrary type T, in general expect float.
        """
        self._name : str = name
        self._val  : T   = val

    @property
    def name(self) -> str: 
        return self._name
    # @property
    # def val(self) -> T: return self._val

    def __repr__(self):
        return f"<Parameter {self.name}={self.get_val()}, at {hex(id(self))}>"

    def get_val(self) -> T:
        return self._val

    def set_val(self, val: T) -> None:
        self._val = val

    

class ParameterManager:
    def __init__(self, 
                 physdata: dict[str, float] | None = None, 
                 scale: float = 4.8):
        """Parameter manager

        Parameters
        ----------
        physdata : dict[str, float], optional
            Dictionary of defualt values (lifetimes, masses), by default DEFAULT_PARAMS which
            is obtained from `flavio.default_parameters.get_central_all()`
        scale : float, optional
            Renorm scale, by default 4.8
        """
        self._params  : dict[str, Parameter]      = {}
        self._scale   : float                     = scale
        self._wc_obj  : flavio.WilsonCoefficients = flavio.WilsonCoefficients()
        self._aliases : dict[str, str]            = {}
        self._load_physdata(physdata if physdata is not None else DEFAULT_PARAMS)

    @property
    def params(self) -> dict[str, Parameter]: return self._params
    @property
    def scale(self) -> float: return self._scale
    @property
    def aliases(self) -> dict[str, str]: return self._aliases

    def __getitem__(self, parname: str) -> Any:
        return self.get_val(parname)

    def get(self, parname: str, default: Any = None) -> Any:
        name = self.aliases.get(parname, parname)
        par = self.params.get(name)
        if par is None:
            return default
        return par.get_val()

    def _load_physdata(self, physdata: dict[str, float]) -> None:
        for iname, ival in physdata.items():
            self.register_param(iname, ival)

    def _update(self, pars: dict[str, Parameter]):
        self.params.update(pars)

    def merge_into(self, other: ParameterManager) -> ParameterManager:
        """Merges two managers, self into other, meaning parameters/aliases of other are
        preserved if they are in both.

        Parameters
        ----------
        other : ParameterManager
            Manager to merge into

        Returns
        -------
        ParameterManager
            Merged parameter manager
        """
        unique_pars = {iname : ipar for iname, ipar in self.params.items() if iname not in other.params}
        unique_aliases = {ialias : iparname for ialias, iparname in self.aliases if ialias not in other.aliases}
        other._update(unique_pars)
        other.aliases.update(unique_aliases)
        return other

    def merge_onto(self, other: ParameterManager) -> ParameterManager:
        """Merges two managers, self onto other, meaning parameters/aliases of self are
        preserved if they are in both.

        Parameters
        ----------
        other : ParameterManager
            Manager to merge onto

        Returns
        -------
        ParameterManager
            Merged parameter manager
        """
        unique_pars = {iname : ipar for iname, ipar in other.params.items() if iname not in self.params}
        unique_aliases = {ialias : iparname for ialias, iparname in other.aliases if ialias not in self.aliases}
        self._update(unique_pars)
        self.aliases.update(unique_aliases)
        return self
    
    # def merge(self, other_manager: ParameterManager) -> None:
    #     if other_manager is self:
    #         return
    #     self.update(other_manager.params)

    def add_param(self, par: Parameter, override = False) -> None:
        if override:
            self._update({par.name : par})
            return
        
        if par.name in self.params:
            dpar = self.get_param(par.name)
            if par is dpar:
                return
            else:
                logger.debug(f"Parameter {par.name} already in manager. This assignment is being skipped.")
            return

        self._update({par.name : par})

    def register_param(self, name: str, val: float) -> None:
        par = Parameter(name, val)
        self.add_param(par)

    def add_alias(self, parname: str, alias: str) -> None:
        if parname not in self.params:
            logger.warning(f"Adding alias for parameter {parname} which is not in manager")
        if alias in self.aliases:
            logger.warning(f"Alias {alias} already exists and is linked to {self.aliases[alias]}. Will be overriden.")
        self.aliases[alias] = parname
    
    def get_val(self, parname: str) -> float:
        return self.get_param(parname).get_val()
    
    def set_val(self, parname: str, val: float) -> None:
        par = self.get_param(parname)
        par.set_val(val)

    def set_vals(self, parvals: dict) -> None:
        for ipar in parvals:
            self.set_val(ipar, parvals[ipar])
    
    def get_param(self, parname: str) -> Parameter:
        name = self.aliases.get(parname, parname)
        return self.params[name]

    def print(self) -> None:
        print(f"ParameterManager at {hex(id(self))}, parameters:")
        for ipar in self.params.values():
            print(f"  {ipar.name:<36} : {ipar}")

    def print_values(self) -> None:
            print(f"ParameterManager at {hex(id(self))}, parameters:")
            for ipar in self.params.values():
                print(f"  {ipar.name:<36} : {ipar.get_val()}")
    
    def set_wc(self, wc: dict[str, complex], eft="WET", basis="flavio") -> None:
        for iwc, ival in wc.items():
            if "WCRe:"+iwc not in self.params:
                self.register_param("WCRe:"+iwc, ival.real)
            else:
                self.set_val("WCRe:"+iwc, ival.real)
            
            if "WCIm:"+iwc not in self.params:
                self.register_param("WCIm:"+iwc, ival.imag)
            else:
                self.set_val("WCIm:"+iwc, ival.imag)
        # Need to use _getWC_values() to ensure still use any values that are unchanged
        self._wc_obj.set_initial(self._get_wc_values(), self.scale, eft, basis)

    def _get_wc_values(self) -> dict[str, complex]:
        wc = {}
        for ikey, ipar in self.params.items():
            if not (ikey[:5] == "WCRe:" or ikey[5:] == "WCIm:"):
                continue
            # Handle Wilson Coefficients
            name_wc = ikey[5:]
            iwc = ipar.get_val() if "WCRe:" in ikey else 1.0j*ipar.get_val()
            if name_wc not in wc:
                wc[name_wc] = 0.0
            wc[name_wc] += iwc
        return wc

    def get_wc(self, eft="WET", basis="flavio") -> flavio.WilsonCoefficients:
        wc = self._get_wc_values()
        self._wc_obj.set_initial(wc, self.scale, eft, basis)
        return self._wc_obj



class ParameterUser:
    _name: str = ""
    def __init__(self):
        self._pm = ParameterManager()
        self._initialize_params(self.pm)

    @property
    def name(self) -> str: return self._name
    @property
    def pm(self) -> ParameterManager: return self._pm

    def define_userparams(self) -> dict[str, Any]:
        # To implement in derived class!
        return {}

    def user_params_defaults(self) -> dict[str, Any]:
        return {f"{self.name}:{ipar}" : ival for ipar, ival in self.define_userparams().items()}

    def _initialize_params(self, manager: ParameterManager) -> None:
        params = copy.deepcopy(self.user_params_defaults())
        for iname, ival in params.items():
            manager.register_param(iname, ival)

    def set_parammanager(self, manager: ParameterManager) -> None:
        self._pm = manager
        self._initialize_params(self.pm)

    # def get_param(self, name: str, prefix: str = None) -> Any:
    #     parname = name if prefix is None else f"{prefix}:{name}"
    #     return self.pm.get_param(parname)

    def get_param(self, name: str) -> Any:
        return self.pm.get_val(name)

    def set_param(self, name: str, val: Any) -> None:
        self.pm.set_val(name, val)

    def get_userparam(self, name: str) -> Any:
        return self.get_param(f"{self.name}:{name}")

    def set_userparam(self, name: str, val: Any) -> None:
        self.pm.set_val(f"{self.name}:{name}", val)