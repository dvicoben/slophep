from math import sqrt
from slophep.Predictions.FormFactorBase import FormFactor

class FormFactorOneHalfpToOneHalfp(FormFactor):
    def __init__(self, 
                 B: str,
                 C: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._C = C
        self._name: str            = "OneHalfpToOneHalfp_FFBase"
        self._internalparams: dict = {
            "Mb"         : self.par[f'm_{self.B}'],
            "Mc"         : self.par[f'm_{self.C}'],
        }

    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The final state hadron"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        return (mB**2 + mV**2 - q2) / (2 * mB * mV)

    def get_ff(self, q2: float) -> dict:
        """Calculate form factors at particular q2.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs
        """
        ff = {
            "Vt"  : 0.0,
            "V0"  : 0.0,
            "Vp"  : 0.0,
            "At"  : 0.0,
            "A0"  : 0.0,
            "Ap"  : 0.0,
            "T0"  : 0.0,
            "T50" : 0.0,
            "Tp"  : 0.0,
            "T5p" : 0.0
        }
        return ff