import flavio



class FormFactor:
    def __init__(self, 
                 par: dict = None, 
                 scale: float = None):
        self._par = flavio.default_parameters.get_central_all() if type(par) == type(None) else par
        self._scale = 4.8 if not scale else scale
        
        self._name = "FFBase"
        self._ffpar = {}
        self._params = []

    @property
    def name(self) -> str: return self._name
    @property
    def scale(self) -> float: return self._scale
    @property
    def par(self) -> dict: return self._par
    @property
    def ffpar(self) -> dict: return self._ffpar
    @property
    def params(self) -> list[str]: return self._params

    def set_ff(self, **kwargs):
        for elem in kwargs:
            if elem not in self.ffpar:
                raise KeyError(f"{elem} is not a parameter for {self.name}")
            if kwargs[elem] == None:
                continue
            self._ffpar[elem] = kwargs[elem]

    def set_ff_fromlist(self, params: list[float]):
        pars = {self.params[k] : params[k] for k in range(len(self.params))}
        self.set_ff(**pars)

    def get_ff(self, q2: float) -> dict:
        return {
            "V"   : 0.0,
            "A1"  : 0.0,
            "A0"  : 0.0,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "A12" : 0.0,
            "T23" : 0.0
        }