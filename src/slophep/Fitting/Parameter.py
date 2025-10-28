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
    def name(self) -> str: return self._name
    @property
    def val(self) -> float: return self._val
    @property
    def limlo(self) -> float: return self._limlo
    @property
    def limhi(self) -> float: return self._limhi
    @property
    def fixed(self) -> bool: return self._fixed
    @property
    def constrain(self) -> bool: return self._constrain
    @property
    def constraintvals(self) -> list[float]: return self._constraintvals
    @property
    def uuid(self) -> str: return self._uuid

    def fix(self, val: float = None):
        if val:
            self._val = val
        self._fixed = True
    
    def setlim(self, lo: float, hi: float):
        if lo:
            self._limlo = lo
        if hi:
            self._limhi = hi

    def constrainto(self, mu: float, sigma: float):
        self._constrain = True
        self._constraintvals = [mu, sigma]
    
    def setVal(self, val: float):
        self._val = val

    def getVal(self):
        return self.val





class ParameterManager:
    def __init__(self):
        self._params = {}
        self._id_map = {}

    @property
    def params(self) -> dict[str, Parameter]: return self._params
    @property
    def id_map(self) -> dict[str, str]: return self._id_map

    def addParam(self, par: Parameter):
        if par.uuid in self.id_map:
            raise KeyError(f"Please change parameter name, Parmeter {par.name} with ID {par.uuid} already exists")
        self._params[par.uuid] = par
        self._id_map[par.name] = par.uuid
    
    def getVal(self, parname: str):
        uuid = self.id_map[parname]
        return self.params[uuid].getVal()
    
    def setVal(self, parname: str, val: float):
        uuid = self.id_map[parname]
        self.params[uuid].setVal(val)

    def setVals(self, parvals: dict):
        for ipar in parvals:
            self.setVal(ipar, parvals[ipar])
    
    def getParam(self, parname: str):
        uuid = self._id_map[parname]
        return self.params[uuid]