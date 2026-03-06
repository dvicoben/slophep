import flavio
from math import sqrt
from slophep.Predictions.FormFactorBase import FormFactor


class FormFactorBToDz(FormFactor):
    def __init__(self, 
                 B: str,
                 C: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._C = C
        self._name: str            = "BToDz_FFBase"
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
        """The Charmed Vector meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mC = self.internalparams["Mc"]
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    def get_ff(self, q2: float) -> dict:
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


class FormFactorBToD1st(FormFactor):
    def __init__(self, 
                 B: str,
                 C: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._C = C
        self._name: str            = "BToD1st_FFBase"
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
        """The Charmed Vector meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mC = self.internalparams["Mc"]
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    def get_ff(self, q2: float) -> dict:
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



class FormFactorBToD1(FormFactor):
    def __init__(self, 
                 B: str,
                 C: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._C = C
        self._name: str            = "BToD1st_FFBase"
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
        """The Charmed Vector meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mC = self.internalparams["Mc"]
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    def get_ff(self, q2: float) -> dict:
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


class FormFactorBToD2st(FormFactor):
    def __init__(self, 
                 B: str,
                 C: str,
                 par: dict = None, 
                 scale: float = None):
        super().__init__(par, scale)
        
        self._B = B
        self._C = C
        self._name: str            = "BToD1st_FFBase"
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
        """The Charmed Vector meson"""
        return self._C
    
    def w(self, q2: float) -> float:
        mB = self.internalparams["Mb"]
        mC = self.internalparams["Mc"]
        return (mB**2 + mC**2 - q2) / (2 * mB * mC)

    def get_ff(self, q2: float) -> dict:
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