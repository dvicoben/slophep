import slophep.Predictions.FormFactorsBToV as FFBToV
import numpy as np

class ISGW2(FFBToV.ISGW2_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_ISGW2"


class BLPR(FFBToV.BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPR"


class BLPR_Hammer(FFBToV.BLPR_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPR_Hammer"


class BLPRXP(FFBToV.BLPRXP_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPRXP"


class CLN(FFBToV.CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN"


class CLN_Hammer(FFBToV.CLN_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN_Hammer"


class CLN2(FFBToV.CLN2_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN2"
        self.internalparams.update({
            "qiqj" : "bc"
        })


class BGL(FFBToV.BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BGL"


class BGL_Hammer(FFBToV.BGL_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BGL_Hammer"


class BSZ(FFBToV.BSZ_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BSZ"


class HPQCD(FFBToV.HPQCD_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_HPQCD"


class BGL_FLabMILC(FFBToV.BGLGeneric_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, 2, 2, 2, 2)
        # Meant to reproduce https://arxiv.org/pdf/2105.14019
        # NOTE: c0 is fixed by kinematical constraint in Eq. (72)
        self._ffpar = {
            "a0": 0.03303572420791651,
            "a1": -0.15617843109815369,
            "a2": -0.12035867559458847,
            "b0": 0.012291940805423925,
            "b1": -0.0033630959482196035,
            "b2": 0.06700760806488626,
            "c0": 0.002058744153712929,
            "c1": -0.005812936725638205,
            "c2": -0.013233078304864712,
            "d0": 0.05086532458850226,
            "d1": -0.3282753381548069,
            "d2": -0.023361655819103194
        }

        internalparams = {
            "chip"       : 5.131e-4,
            "chim"       : 3.894e-4,
            "chimL"      : 1.9421e-2,
            "BcStatesf"  : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStatesg"  : np.array([6.329, 6.920, 7.020]),
            "BcStatesP1" : np.array([6.275, 6.842, 7.250])
        }
        self.internalparams.update(internalparams)


class BGL_JLQCD(FFBToV.BGLGeneric_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, 2, 2, 2, 2)
        # Meant to reproduce https://arxiv.org/pdf/2306.05657
        # NOTE: c0 is fixed by kinematical constraint in Eq. (33)
        self._ffpar = {
            "a0" : 0.0291,
            "a1" : -0.045,
            "a2" : -1.0,
            "b0" : 0.01198,
            "b1" : 0.018,
            "b2" : -0.10,
            "c0" : 0.002006,
            "c1" : 0.0013,
            "c2" : -0.03,
            "d0" : 0.0484,
            "d1" : -0.059,
            "d2" : -0.9
        }

        internalparams = {
            "chip"       : 5.131e-4,
            "chim"       : 3.894e-4,
            "chimL"      : 1.9421e-2,
            "BcStatesf"  : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStatesg"  : np.array([6.329, 6.920, 7.020]),
            "BcStatesP1" : np.array([6.275, 6.842, 7.250])
        }
        self.internalparams.update(internalparams)