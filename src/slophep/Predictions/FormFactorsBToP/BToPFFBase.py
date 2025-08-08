import flavio
from math import sqrt
from slophep.Predictions.FormFactorBase import FormFactor

class FormFactorBToP(FormFactor):
    def __init__(self, 
                 B: str,
                 P: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._P = P
        self._name: str            = "BToPFFBase"
        self._internalparams: dict = {
            "Mb"         : self.par[f'm_{self.B}'],
            "Mc"         : self.par[f'm_{self.P}'],
        }

    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def P(self) -> str:
        """The Charmed Vector meson"""
        return self._P
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        return (mB**2 + mV**2 - q2) / (2 * mB * mV)

    def get_ff(self, q2: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.
        Must return in basis f+, f0, fT

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs 
        """
        return {
            "f+" : 0.0,
            "f0" : 0.0,
            "fT" : 0.0
        }
    