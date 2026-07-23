from typing import Any
import numpy as np

import slophep.FormFactors.FormFactorsBToV as FFBToV
from slophep.Tools.SamplingTools import fluctsettings, FluctType


class ISGW2(FFBToV.FFBToV_ISGW2):
    _name = "FFBdToDst@ISGW2"
    def __init__(self):
        super().__init__("B0", "D*+")


class BLPR(FFBToV.FFBToV_BLPR):
    _name = "FFBdToDst@BLPR"
    def __init__(self):
        super().__init__("B0", "D*+")


class BLPRXP(FFBToV.FFBToV_BLPRXP):
    _name = "FFBdToDst@BLPRXP"
    def __init__(self):
        super().__init__("B0", "D*+")
    

class CLN(FFBToV.FFBToV_CLN):
    _name = "FFBdToDst@CLN"
    def __init__(self):
        super().__init__("B0", "D*+")


class CLN2(FFBToV.FFBToV_CLN2):
    _name = "FFBdToDst@CLN2"
    def __init__(self):
        super().__init__("B0", "D*+")


class BGL(FFBToV.FFBToV_BGL):
    _name = "FFBdToDst@BGL"
    def __init__(self):
        super().__init__("B0", "D*+")


class BSZ(FFBToV.FFBToV_BSZ):
    _name = "FFBdToDst@BSZ"
    def __init__(self):
        super().__init__("B0", "D*+")


class HPQCD(FFBToV.FFBToV_HPQCD):
    _name = "FFBdToDst@HPQCD2023"
    def __init__(self):
        super().__init__("B0", "D*+")


class BGL_FNALMILC(FFBToV.FFBToV_BGLGeneric):
    _name = "FFBdToDst@BGL_FNAL/MILC2021"
    def __init__(self):
        super().__init__("B0", "D*+", 3, 3, 3, 3, 0, 0, 0)

    def define_userparams(self) -> dict[str, Any]:
        # Meant to reproduce https://arxiv.org/pdf/2105.14019
        # NOTE: c0 is fixed by kinematical constraint in Eq. (72)
        pars = super().define_userparams()
        ffpar = {
            "a_g_0"      : 0.03303572420791651   ,
            "a_g_1"      : -0.15617843109815369  ,
            "a_g_2"      : -0.12035867559458847  ,
            "a_f_0"      : 0.012291940805423925  ,
            "a_f_1"      : -0.0033630959482196035,
            "a_f_2"      : 0.06700760806488626   ,
            "a_F1_0"     : 0.002058744153712929  ,
            "a_F1_1"     : -0.005812936725638205 ,
            "a_F1_2"     : -0.013233078304864712 ,
            "a_F2_0"     : 0.05086532458850226   ,
            "a_F2_1"     : -0.3282753381548069   ,
            "a_F2_2"     : -0.023361655819103194 ,
            "chig"       : 5.131e-4 ,
            "chif"       : 3.894e-4 ,
            "chiF1"      : 3.894e-4 ,
            "chiF2"      : 1.9421e-2,
            "BcStates1p" : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStates1m" : np.array([6.329, 6.920, 7.020]),
            "BcStates0m" : np.array([6.275, 6.842, 7.250])
        }
        pars.update(ffpar)
        return pars



class BGL_JLQCD(FFBToV.FFBToV_BGLGeneric):
    _name = "FFBdToDst@BGL_JLQCD2023"
    def __init__(self):
        super().__init__("B0", "D*+", 3, 3, 3, 3, 0, 0, 0)

    def define_userparams(self) -> dict[str, Any]:
        # Meant to reproduce https://arxiv.org/pdf/2105.14019
        # NOTE: c0 is fixed by kinematical constraint in Eq. (72)
        pars = super().define_userparams()
        ffpar = {
            "a_g_0"      : 0.0291   ,
            "a_g_1"      : -0.045   ,
            "a_g_2"      : -1.0     ,
            "a_f_0"      : 0.01198  ,
            "a_f_1"      : 0.018    ,
            "a_f_2"      : -0.10    ,
            "a_F1_0"     : 0.002006 ,
            "a_F1_1"     : 0.0013   ,
            "a_F1_2"     : -0.03    ,
            "a_F2_0"     : 0.0484   ,
            "a_F2_1"     : -0.059   ,
            "a_F2_2"     : -0.9     ,
            "chig"       : 5.131e-4 ,
            "chif"       : 3.894e-4 ,
            "chiF1"      : 3.894e-4 ,
            "chiF2"      : 1.9421e-2,
            "BcStates1p" : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStates1m" : np.array([6.329, 6.920, 7.020]),
            "BcStates0m" : np.array([6.275, 6.842, 7.250])
        }
        pars.update(ffpar)
        return pars



class BGL_Hammer(FFBToV.FFBToV_BGLGeneric):
    _name = "FFBdToDst@BGL_Hammer"
    def __init__(self):
        super().__init__("B0", "D*+", 3, 3, 3, 2, 0, 0, 0)

    def define_userparams(self) -> dict[str, Any]:
        # Meant to reproduce https://arxiv.org/pdf/2105.14019
        # NOTE: c0 is fixed by kinematical constraint in Eq. (72)
        pars = super().define_userparams()
        ffpar = {
            "a_g_0"  : 0.00038   ,
            "a_g_1"  : 0.026905  ,
            "a_g_2"  : 0.        ,
            "a_f_0"  : 0.00055   ,
            "a_f_1"  : -0.0020370,
            "a_f_2"  : 0.        ,
            "a_F1_1" : -0.000433 ,
            "a_F1_2" : 0.005353  ,
            "a_F2_0" : 0.007     ,
            "a_F2_1" : -0.036    ,
            # internalparams
            "Vcb"        : 41.5e-3 ,                       
            "chiF1"      : 3.068e-4, # GeV^-2
            "chif"       : 3.068e-4, # GeV^-2
            "chig"       : 5.280e-4, # GeV^-2
            "chiF2"      : 2.466e-3,
            "nc"         : 2.6     ,
            "etaEW"      : 1.0066  ,
            "BcStates1p" : np.array([6.730, 6.736, 7.135, 7.142]),  # GeV
            "BcStates1m" : np.array([6.337, 6.899, 7.012, 7.280]),  # GeV
            "BcStates0m" : np.array([6.275, 6.842, 7.250])          # GeV
        }
        pars.update(ffpar)
        return pars

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict:
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags
        
        Note that there is an additional 1./(etaEW*Vcb) factor applied to FFs as in Hammer.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        etaEWVcb = self.get_userparam("etaEW")*self.get_userparam("Vcb")
        ff = super().get_ff(q2)
        return {k : ff[k]/etaEWVcb for k in ff}


class BLPR_Hammer(FFBToV.FFBToV_BLPR):
    _name = "FFBToDst@BLPR_Hammer"
    def __init__(self):
        super().__init__("B0", "D*+")

    def define_userparams(self):
        ffpar = {
            "RhoSq"     : 1.24      ,
            "Chi21"     : -0.06     ,
            "Chi2p"     : 0.0       ,
            "Chi3p"     : 0.05      ,
            "Eta1"      : 0.30      ,
            "Etap"      : -0.05     ,
            "dV20"      : 0.0       ,
            # internalparams:
            "mb"        : 4.710     ,
            "delta_mbc" : 3.4       ,
            "normscale" : 1.0       ,
            "ash"       : 0.26/np.pi,
            "ebReb"     : 0.861     ,
            "ecRec"     : 0.822     ,
            "mbarB"     : 5.313     ,
            "lam1"      : -0.3      , # GeV^2
        }
        return ffpar


class CLN_Hammer(FFBToV.FFBToV_CLN):
    _name = "FFBdToDst@CLN_Hammer"
    def __init__(self):
        super().__init__("B0", "D*+")

    def define_userparams(self):
        ffpar = {
            "RhoSq" : 1.207,
            "h_A1"  : 0.908,
            "R1"    : 1.401,
            "R2"    : 0.854,
            "R0"    : 1.15
        }
        return ffpar