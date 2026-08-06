"""
Bs->Ds1 Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD1_ISGW2):
    _name = "FFBsToDs1@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds1+")

    def define_userparams(self):
        ffpar = {
            "msb"     : 5.2                    ,
            "msd"     : 0.55                   ,
            "bb2"     : 0.54*0.54              ,
            "mbb"     : 5.38                   ,
            "msq"     : 1.82                   ,
            "bx2"     : 0.41*0.41              ,
            "mbx"     : (5.0*2.61+3.0*2.54)/8.0,
            "mqm"     : 0.1                    ,
            "nfp"     : 3.0                    ,
            "SmearQ2" : True
        }
        return ffpar


@FFregistry.register
class LLSW(FFBToDstst.FFBToD1_LLSW):
    _name = "FFBsToDs1@LLSW"
    def __init__(self):
        super().__init__("Bs", "Ds1+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD1_BLR):
    _name = "FFBsToDs1@BLR"
    def __init__(self):
        super().__init__("Bs", "Ds1+")