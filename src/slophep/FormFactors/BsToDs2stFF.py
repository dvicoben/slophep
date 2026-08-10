"""
Bs->Ds2* Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD2st_ISGW2):
    _name = "FFBsToDs2st@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds2*+")

    def define_userparams(self):
        ffpar = {
            "msb" : 5.2                    ,
            "msd" : 0.55                   ,
            "bb2" : 0.54*0.54              ,
            "mbb" : 5.38                   ,
            "msq" : 1.82                   ,
            "bx2" : 0.41*0.41              ,
            "mbx" : (5.0*2.61+3.0*2.54)/8.0,
            "mqm" : 0.1                    ,
            "nfp" : 3.0
        }
        return ffpar


@FFregistry.register
class LLSW(FFBToDstst.FFBToD2st_LLSW):
    _name = "FFBsToDs2st@LLSW"
    def __init__(self):
        super().__init__("Bs", "Ds2*+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD2st_BLR):
    _name = "FFBsToDs2st@BLR"
    def __init__(self):
        super().__init__("Bs", "Ds2*+")