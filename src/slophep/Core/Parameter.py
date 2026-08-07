from __future__ import annotations
import flavio
from typing import TypeVar, Any

from slophep.Core.physdata import DEFAULT_PARAMS
import logging

logger = logging.getLogger(__name__)

T = TypeVar("T")

class Parameter:
    def __init__(self, 
                 name: str,
                 val: T):
        """Parameter object

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
        """Parameter name"""
        return self._name
    # @property
    # def val(self) -> T: return self._val

    def __repr__(self):
        return f"<Parameter {self.name}={self.get_val()}, at {hex(id(self))}>"

    def get_val(self) -> T:
        """Get the parameter value"""
        return self._val

    def set_val(self, val: T) -> None:
        """Set the parameter value"""
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
        self._wc_obj  : flavio.WilsonCoefficients = flavio.WilsonCoefficients()
        self._aliases : dict[str, str]            = {}
        # Load default params - masses, lifetimes, etc.
        self._load_physdata(physdata if physdata is not None else DEFAULT_PARAMS)
        # Add scale parameter
        self.register_param("SCALE", scale)

    @property
    def params(self) -> dict[str, Parameter]: 
        """Dictionary of parameters"""
        return self._params
    @property
    def scale(self) -> float: 
        """Renorm scale"""
        return self.get_val("SCALE")
    @property
    def aliases(self) -> dict[str, str]: 
        """Dictionary of aliases form {alias : qualified_param_name}"""
        return self._aliases

    def set_scale(self, scale: float) -> None:
        """Set the scale

        Parameters
        ----------
        scale : float
        """
        self.set_val("SCALE", scale)

    def __getitem__(self, parname: str) -> Any:
        return self.get_val(parname)

    def __contains__(self, parname: str) -> bool:
        return (parname in self.params) or (parname in self.aliases)

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
        preserved if they are in both. Note this alters and returns other, it does not create a new manager.

        Parameters
        ----------
        other : ParameterManager
            Manager to merge into

        Returns
        -------
        ParameterManager
            Merged parameter manager, other with the contents of self without overwriting other
        """
        unique_pars = {iname : ipar for iname, ipar in self.params.items() if iname not in other.params}
        unique_aliases = {ialias : iparname for ialias, iparname in self.aliases if ialias not in other.aliases}
        other._update(unique_pars)
        other.aliases.update(unique_aliases)
        return other

    def merge_onto(self, other: ParameterManager) -> ParameterManager:
        """Merges two managers, self onto other, meaning parameters/aliases of self are
        preserved if they are in both. Note this alters and returns self, it does not create a new manager.

        Parameters
        ----------
        other : ParameterManager
            Manager to merge onto

        Returns
        -------
        ParameterManager
            Merged parameter manager, self with the contents of other without overwriting self
        """
        unique_pars = {iname : ipar for iname, ipar in other.params.items() if iname not in self.params}
        unique_aliases = {ialias : iparname for ialias, iparname in other.aliases if ialias not in self.aliases}
        self._update(unique_pars)
        self.aliases.update(unique_aliases)
        return self

    def add_param(self, par: Parameter, override = False) -> None:
        """Add parameter to manager. If it exists already it will not be added.

        Parameters
        ----------
        par : Parameter
            Parameter to add
        override : bool, optional
            Whether to ignore duplication checks and override existing parameter, by default False
        """
        if override:
            if par.name in self.aliases:
                logger.warning(f"Parameter {par.name} registered as an alias {par.name}->{self.aliases[par.name]}. Alias will be removed.")
                self.remove_alias(par.name)
            self._update({par.name : par})
            return
        
        if par.name in self.aliases:
            logger.warning(f"Parameter {par.name} registered as an alias {par.name}->{self.aliases[par.name]}. This assignment is being skipped.")
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
        """Create a parameter and add it to manager

        Parameters
        ----------
        name : str
            Name of parameter
        val : float
            Value of parameter
        """
        par = Parameter(name, val)
        self.add_param(par)

    def add_alias(self, parname: str, alias: str) -> None:
        """Add an alias to a parameter, can use in case want several names to refer to the same parameter

        Parameters
        ----------
        parname : str
            Qualified parameter name to create alias for
        alias : str
            Alias for the parameter
        """
        if parname not in self.params:
            logger.warning(f"Adding alias for parameter {parname} which is not in manager")
        if alias in self.aliases:
            logger.warning(f"Alias {alias} already exists and is linked to {self.aliases[alias]}. Will be overriden.")
        if alias in self.params:
            logger.error(f"Cannot add alias for an existing parameter, alias not added")
            return
        self.aliases[alias] = parname

    def remove_alias(self, alias: str) -> None:
        """Remove an alias

        Parameters
        ----------
        alias : str
            Alias to remove
        """
        parname = self.aliases.pop(alias)
        logger.info(f"Removed alias {alias} linked to {parname}")
    
    def get_val(self, parname: str) -> Any:
        """Get value for a particular parameter

        Parameters
        ----------
        parname : str
            Parameter name or alias

        Returns
        -------
        Any
            Value of the parameter
        """
        return self.get_param(parname).get_val()
    
    def set_val(self, parname: str, val: Any) -> None:
        """Set value of a parameter

        Parameters
        ----------
        parname : str
            Parameter name or alias
        val : Any
            Value to set
        """
        par = self.get_param(parname)
        par.set_val(val)

    def set_vals(self, parvals: dict) -> None:
        """Set values for series of parameters

        Parameters
        ----------
        parvals : dict
            Parameters and values to set them to.
        """
        for ipar in parvals:
            self.set_val(ipar, parvals[ipar])
    
    def get_param(self, parname: str) -> Parameter:
        """Get a parameter object

        Parameters
        ----------
        parname : str
            Name or alias of parameter

        Returns
        -------
        Parameter
            _description_
        """
        name = self.aliases.get(parname, parname)
        return self.params[name]

    def print(self) -> None:
        """Print manager contents
        """
        print(f"ParameterManager at {hex(id(self))}, parameters:")
        for ipar in self.params.values():
            print(f"  {ipar.name:<36} : {ipar}")
        for ialias, ipname in self.aliases.items():
            print(f"  {ialias:<36} : ALIAS->{ipname}")

    def print_values(self) -> None:
        """Print values of parameters in manager
        """
        print(f"ParameterManager at {hex(id(self))}, parameters:")
        for ipar in self.params.values():
            print(f"  {ipar.name:<36} : {ipar.get_val()}")
        for ialias, ipname in self.aliases.items():
            print(f"  {ialias:<36} : ALIAS->{ipname}")
    
    def set_wc(self, wc: dict[str, complex], eft: str = "WET", basis: str = "flavio") -> None:
        """Set wilson coefficient values

        Parameters
        ----------
        wc : dict[str, complex]
            Wilson coefficients, in flavio basis
        eft : str, optional
            EFT, by default "WET"
        basis : str, optional
            WC vasis, by default "flavio"
        """
        for iwc, ival in wc.items():
            if "WCRe:"+iwc not in self.params:
                self.register_param("WCRe:"+iwc, ival.real)
            else:
                self.set_val("WCRe:"+iwc, ival.real)
            
            if "WCIm:"+iwc not in self.params:
                self.register_param("WCIm:"+iwc, ival.imag)
            else:
                self.set_val("WCIm:"+iwc, ival.imag)
        # Need to use _get_wc_values() to ensure still use any values that are unchanged
        self._wc_obj.set_initial(self._get_wc_values(), self.scale, eft, basis)

    def _get_wc_values(self) -> dict[str, complex]:
        wc = {}
        for ikey, ipar in self.params.items():
            if not (ikey[:5] == "WCRe:" or ikey[:5] == "WCIm:"):
                continue
            # Handle Wilson Coefficients
            name_wc = ikey[5:]
            iwc = ipar.get_val() if "WCRe:" in ikey else 1.0j*ipar.get_val()
            if name_wc not in wc:
                wc[name_wc] = 0.0
            wc[name_wc] += iwc
        return wc

    def get_wc(self, eft: str ="WET", basis: str = "flavio") -> flavio.WilsonCoefficients:
        """Get wilson coefficient object

        Parameters
        ----------
        eft : str, optional
            EFT, by default "WET"
        basis : str, optional
            WC basis, by default "flavio"

        Returns
        -------
        flavio.WilsonCoefficients
            WC object
        """
        wc = self._get_wc_values()
        self._wc_obj.set_initial(wc, self.scale, eft, basis)
        return self._wc_obj



class ParameterUser:
    _name: str = ""
    def __init__(self):
        """Base clase for objects needing access to a ParameterManager
        """
        self._pm = ParameterManager()
        self._initialize_params(self.pm)

    @property
    def name(self) -> str: 
        """User name, should be class specific, not instance specific"""
        return self._name
    @property
    def pm(self) -> ParameterManager: 
        """User parameter manager"""
        return self._pm

    def define_userparams(self) -> dict[str, Any]:
        """Dictionary of parameters relevant for this specific user. To implement in child classes.
        
        This method return a dictionary of the form {"paramname" : value}.
        On initialisation the user will register these parameters as {"username:paramname" : value} into the manager.

        Returns
        -------
        dict[str, Any]
            Dictionary of the form {"paramname" : value}
        """
        # To implement in derived class!
        return {}

    def user_params_defaults(self) -> dict[str, Any]:
        """Returns dictionary of user defaults with fully qualified names

        Returns
        -------
        dict[str, Any]
            Dictionary of the form {"username:paramname" : value}
        """
        return {f"{self.name}:{ipar}" : ival for ipar, ival in self.define_userparams().items()}

    def _initialize_params(self, manager: ParameterManager) -> None:
        params = self.user_params_defaults()
        for iname, ival in params.items():
            manager.register_param(iname, ival)

    def set_parammanager(self, manager: ParameterManager) -> None:
        """Set user parameter manager. Registers user parameters if they are not in the manager.

        Parameters
        ----------
        manager : ParameterManager
        """
        self._pm = manager
        self._initialize_params(self.pm)

    # def get_param(self, name: str, prefix: str = None) -> Any:
    #     parname = name if prefix is None else f"{prefix}:{name}"
    #     return self.pm.get_param(parname)

    def get_param(self, name: str) -> Any:
        """Get value for parameter using its full qualified name

        Parameters
        ----------
        name : str
            Parameter qualified name

        Returns
        -------
        Any
            parameter value
        """
        return self.pm.get_val(name)

    def set_param(self, name: str, val: Any) -> None:
        """Set parameter using its full qualified name

        Parameters
        ----------
        name : str
            Parameter qualified name
        val : Any
            Value to set
        """
        self.pm.set_val(name, val)

    def get_userparam(self, name: str) -> Any:
        """Get value user parameter using its internal name, 
        i.e. provided `name`, searches for `USERNAME:name`

        Parameters
        ----------
        name : str
            Parameter "unqualified" name

        Returns
        -------
        Any
            Parameter value
        """
        return self.get_param(f"{self.name}:{name}")

    def set_userparam(self, name: str, val: Any) -> None:
        """Set user parameter using its internal name, 
        i.e. provided `name`, sets for `USERNAME:name`

        Parameters
        ----------
        name : str
            Parameter "unqualified" name
        val : Any
            Value to set
        """
        self.pm.set_val(f"{self.name}:{name}", val)