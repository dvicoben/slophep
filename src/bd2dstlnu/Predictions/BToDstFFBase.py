import flavio



class FormFactor:
    def __init__(self, 
                 par: dict = None, 
                 scale: float = None):
        self._par: dict            = flavio.default_parameters.get_central_all() if type(par) == type(None) else par
        self._scale: float         = 4.8 if not scale else scale
        
        self._name: str            = "FFBase"
        self._ffpar: dict          = {}
        self._params: list         = []
        self._internalparams: dict = {}

    @property
    def name(self) -> str: 
        """Name of FF scheme"""
        return self._name
    @property
    def scale(self) -> float:
        """Renorm scale"""
        return self._scale
    @property
    def par(self) -> dict:
        """Dictionary of flavio parameters, defaults to `flavio.default_parameters.get_central_all()`"""
        return self._par
    @property
    def ffpar(self) -> dict: 
        """Form factor (floating) parameters, make sure to overwrite in derived classes"""
        return self._ffpar
    @property
    def params(self) -> list[str]: 
        """List of (floating) form factor parameters (the ones in ffpar), make sure to overwrite in derived classes"""
        return self._params
    @property
    def internalparams(self) -> dict: 
        """Internal parameters used in FF computation (e.g. alpha_s, quark/meson masses), 
        where possible take from FormFactor.par, but things may differ if an implementation is taken from HAMMER
        so this is kept seperate"""
        return self._internalparams

    def set_ff(self, **kwargs):
        """Set form factors - argument names should be as in FormFactor.params"""
        for elem in kwargs:
            if elem not in self.ffpar:
                raise KeyError(f"{elem} is not a parameter for {self.name}")
            if kwargs[elem] == None:
                continue
            self._ffpar[elem] = kwargs[elem]

    def set_ff_fromlist(self, params: list[float]):
        """Set form factors

        Parameters
        ----------
        params : list[float]
            Form factor parameters, in order as in FormFactor.params
        """
        pars = {self.params[k] : params[k] for k in range(len(self.params))}
        self.set_ff(**pars)

    def get_ff(self, q2: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.
        Must return in basis V, A0, A1, A12, T1, T2, T23

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs V, A0, A1, A12, T1, T2, T23
        """
        return {
            "V"   : 0.0,
            "A1"  : 0.0,
            "A0"  : 0.0,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "A12" : 0.0,
            "T23" : 0.0
        }