from slophep.FormFactors.FormFactorBase import FormFactor
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

class FormFactorOneHalfpToOneHalfp(FormFactor):
    _name = "FFBaryonic1/2To1/2@Base"
    def __init__(self, 
                 B: str,
                 C: str):
        """Baryonic 1/2+ -> 1/2+ form-factors"""
        self._B = B
        self._C = C
        super().__init__()
                
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The final state hadron"""
        return self._C
    
    def w(self, q2: float) -> float:
            mB = self.get_param(f"m_{self.B}")
            mP = self.get_param(f"m_{self.P}")
            return (mB**2 + mP**2 - q2) / (2 * mB * mP)
    
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