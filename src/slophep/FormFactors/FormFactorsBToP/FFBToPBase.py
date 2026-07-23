from math import sqrt
from slophep.FormFactors.FormFactorBase import FormFactor
from slophep.Tools.SamplingTools import fluctsettings, FluctType

class FormFactorBToP(FormFactor):
    _name = "FFBToP@Base"
    def __init__(self, B: str, P: str):
        self._B = B
        self._P = P
        super().__init__()
    
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def P(self) -> str:
        """The P meson"""
        return self._P
    
    def w(self, q2: float) -> float:
        mB = self.get_param(f"m_{self.B}")
        mP = self.get_param(f"m_{self.P}")
        return (mB**2 + mP**2 - q2) / (2 * mB * mP)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
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
