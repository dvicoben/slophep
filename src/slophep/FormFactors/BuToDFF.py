"""
B+->D0 Form-factors
"""
import numpy as np
import slophep.FormFactors.FormFactorsBToP as FFBToP
from slophep.Core.user_registry import FFregistry

import logging
logger = logging.getLogger(__name__)

@FFregistry.register
class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBuToD@BSZ"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BLPR(FFBToP.FFBToP_BLPR):
    _name = "FFBuToD@BLPR"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BLPRXP(FFBToP.FFBToP_BLPRXP):
    _name = "FFBuToD@BLPRXP"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BGL(FFBToP.FFBToP_BGL):
    _name = "FFBuToD@BGL"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class CLN(FFBToP.FFBToP_CLN):
    _name = "FFBuToD@CLN"
    def __init__(self):
        super().__init__("B+", "D0")


# Hammer equivalent classes
@FFregistry.register
class CLN_Hammer(FFBToP.FFBToP_CLN):
    _name = "FFBuToD@CLN_Hammer"
    def __init__(self):
        super().__init__("B+", "D0")

    def define_userparams(self):
        ffpar = {
            "RhoSq"  : 1.186,
            "G1"     : 1.082,
            "Delta"  : 1.0  ,
            "a"      : 1.0  ,  # zero recoil expansion
            "mcOnmb" : 0.29 ,
            "lamBar" : 0.48 ,
            "ash"    : 0.26/np.pi
        }
        return ffpar


@FFregistry.register
class BGL_Hammer(FFBToP.FFBToP_BGLGeneric):
    _name = "FFBuToD@BGL_Hammer"
    def __init__(self):
        super().__init__("B+", "D0", 4, 4)

    def define_userparams(self):
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f+_0" : 0.01565,
            "a_f+_1" : -0.0353,
            "a_f+_2" : -0.043,
            "a_f+_3" : 0.194,
            "a_f0_0" : 0.07932,
            "a_f0_1" : -0.214,
            "a_f0_2" : 0.17,
            "a_f0_3" : -0.958,
            #internalparams
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        })
        return ffpar
    

@FFregistry.register
class BLPR_Hammer(FFBToP.FFBToP_BLPR):
    _name = "FFBuToD@BLPR_Hammer"
    def __init__(self):
        logger.info("B->P BLPR_Hammer basis differs in fT from standard SLOP basis, use with care for predictions.")
        super().__init__("B+", "D0")

    def define_userparams(self):
        ffpar = {
            "RhoSq"     : 1.24,
            "Chi21"     : -0.06,
            "Chi2p"     : 0.0,
            "Chi3p"     : 0.05,
            "Eta1"      : 0.30,
            "Etap"      : -0.05,
            "dV20"      : 0.0,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "normscale" : 1.0,
            "ash"       : 0.26/np.pi,
            "mbarB"     : 5.313,
            "lam1"      : -0.3 , # GeV^2
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
        }
        return ffpar