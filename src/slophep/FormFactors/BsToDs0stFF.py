import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD0st_ISGW2):
    _name = "FFBsToDs0st@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds0*+")

    def define_userparams(self):
        ffpar = {
            "msb" : 5.2                ,
            "msd" : 0.55               ,
            "bb2" : 0.54*0.54          ,
            "mbb" : 5.38               ,
            "msq" : 1.82               ,
            "bx2" : 0.41*0.41          ,
            "mbx" : (3.0*2.54+2.46)/4.0,
            "mqm" : 0.1                ,
            "nfp" : 3.0
        }
        return ffpar


@FFregistry.register
class LLSW(FFBToDstst.FFBToD0st_LLSW):
    _name = "FFBsToDs0st@LLSW"
    def __init__(self):
        super().__init__("Bs", "Ds0*+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD0st_BLR):
    _name = "FFBsToDs0st@BLR"
    def __init__(self):
        super().__init__("Bs", "Ds0*+")