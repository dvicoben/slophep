
from math import sqrt
from slophep.FormFactors.FormFactorBase import FormFactor
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

class FormFactorBToD0st(FormFactor):
    _name = "FFBToD0st@Base"
    def __init__(self, B: str, C: str):
        """B->D_0^{*(0/+)} form-factors"""
        self._B = B
        self._C = C
        super().__init__()
    
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The Charmed meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.get_param(f"m_{self.B}")
        mC = self.get_param(f"m_{self.C}")
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.

        Parameters
        ----------
        q2 : float
        """
        return {
            "gP" : 0.0,
            "g+" : 0.0,
            "g-" : 0.0,
            "gT" : 0.0
        }

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
        return self.get_ff_mmeson(q2, self.get_param(f"m_{self.C}"), self.get_param(f"m_{self.B}"))



class FormFactorBToD1st(FormFactor):
    _name = "FFBToD1st@Base"
    def __init__(self, B: str, C: str):
        """B->D1* (i.e. D1') form-factors"""
        self._B = B
        self._C = C
        super().__init__()
    
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The Charmed meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.get_param(f"m_{self.B}")
        mC = self.get_param(f"m_{self.C}")
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.

        Parameters
        ----------
        q2 : float
        """
        return {
            "gS"  : 0.0,
            "gV1" : 0.0,
            "gV2" : 0.0,
            "gV3" : 0.0,
            "gA"  : 0.0,
            "gT1" : 0.0,
            "gT2" : 0.0,
            "gT3" : 0.0
        }

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
        return self.get_ff_mmeson(q2, self.get_param(f"m_{self.C}"), self.get_param(f"m_{self.B}"))



class FormFactorBToD1(FormFactor):
    _name = "FFBToD1@Base"
    def __init__(self, B: str, C: str):
        """B->D1 form-factors"""
        self._B = B
        self._C = C
        super().__init__()
    
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The Charmed meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.get_param(f"m_{self.B}")
        mC = self.get_param(f"m_{self.C}")
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.

        Parameters
        ----------
        q2 : float
        """
        return {
            "fS"  : 0.0,
            "fV1" : 0.0,
            "fV2" : 0.0,
            "fV3" : 0.0,
            "fA"  : 0.0,
            "fT1" : 0.0,
            "fT2" : 0.0,
            "fT3" : 0.0
        }

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
        return self.get_ff_mmeson(q2, self.get_param(f"m_{self.C}"), self.get_param(f"m_{self.B}"))



class FormFactorBToD2st(FormFactor):
    _name = "FFBToD2st@Base"
    def __init__(self, B: str, C: str):
        """B->D2* form-factors"""
        self._B = B
        self._C = C
        super().__init__()
    
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def C(self) -> str:
        """The Charmed meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.get_param(f"m_{self.B}")
        mC = self.get_param(f"m_{self.C}")
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float) -> dict:
        """Calculate form factors at particular q2. To implement in derived class.

        Parameters
        ----------
        q2 : float
        """
        return {
            "kP"  : 0.0,
            "kA1" : 0.0,
            "kA2" : 0.0,
            "kA3" : 0.0,
            "kV"  : 0.0,
            "kT1" : 0.0,
            "kT2" : 0.0,
            "kT3" : 0.0
        }

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
        return self.get_ff_mmeson(q2, self.get_param(f"m_{self.C}"), self.get_param(f"m_{self.B}"))