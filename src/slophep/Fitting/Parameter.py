import uuid

class Parameter:
    def __init__(self, 
                 name: str,
                 val: float, 
                 limlo: float = None, 
                 limhi: float = None,
                 fixed: bool = False):

        self._name = name
        self._val = val
        self._limlo = limlo
        self._limhi = limhi
        self._fixed = fixed
        self._constrain = False
        self._constraintvals = []
        self._uuid = uuid.uuid5(uuid.NAMESPACE_DNS, name).hex

    @property
    def name(self) -> str: 
        """Parameter name"""
        return self._name
    @property
    def val(self) -> float: 
        """Parameter value"""
        return self._val
    @property
    def limlo(self) -> float: 
        """Lower limit"""
        return self._limlo
    @property
    def limhi(self) -> float: 
        """Upper limit"""
        return self._limhi
    @property
    def fixed(self) -> bool: 
        """Whether parameter is fixed"""
        return self._fixed
    @property
    def constrain(self) -> bool: 
        """Whether parameter is constrained"""
        return self._constrain
    @property
    def constraintvals(self) -> list[float]: 
        """The values for gaussian constraint, in form [mu, sigma]"""
        return self._constraintvals
    @property
    def uuid(self) -> str: 
        """Unique hash ID"""
        return self._uuid

    def fix(self, val: float = None):
        """Fix parameter to val"""
        if val:
            self._val = val
        self._fixed = True
    
    def setlim(self, lo: float, hi: float):
        """Set parameter limits"""
        if lo:
            self._limlo = lo
        if hi:
            self._limhi = hi

    def constrainto(self, mu: float, sigma: float):
        """Set parameter constraints"""
        self._constrain = True
        self._constraintvals = [mu, sigma]
    
    def setVal(self, val: float):
        """Set parameter value"""
        self._val = val

    def getVal(self):
        """Return parameter value"""
        return self.val





class ParameterManager:
    def __init__(self):
        self._params = {}
        self._id_map = {}

    @property
    def params(self) -> dict[str, Parameter]: 
        """Dictionary of all parameters in manager, indexed by their uniue hash ID"""
        return self._params
    @property
    def id_map(self) -> dict[str, str]: 
        """Correspondence of parameter names to hash ID"""
        return self._id_map

    def addParam(self, par: Parameter):
        """Add parameter to manager

        Parameters
        ----------
        par : Parameter

        Raises
        ------
        KeyError
            If there is already an existing parameter with the same name
        """
        if par.uuid in self.id_map:
            raise KeyError(f"Please change parameter name, Parmeter {par.name} with ID {par.uuid} already exists")
        self._params[par.uuid] = par
        self._id_map[par.name] = par.uuid
    
    def getVal(self, parname: str):
        """Get value of a parameter

        Parameters
        ----------
        parname : str
            The name of the parameter
        """
        uuid = self.id_map[parname]
        return self.params[uuid].getVal()
    
    def setVal(self, parname: str, val: float):
        """Set value of a parameter

        Parameters
        ----------
        parname : str
            The name of the parameter
        """
        uuid = self.id_map[parname]
        self.params[uuid].setVal(val)

    def setVals(self, parvals: dict):
        """Set value of several parameters

        Parameters
        ----------
        parvals : dict
            Dictionary in form {parameter_name : value}
        """
        for ipar in parvals:
            self.setVal(ipar, parvals[ipar])
    
    def getParam(self, parname: str) -> Parameter:
        """Get Parameter object from manager

        Parameters
        ----------
        parname : str
            Name of the parameter

        Returns
        -------
        Parameter
        """
        uuid = self._id_map[parname]
        return self.params[uuid]