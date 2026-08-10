"""
B+->D2* Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD2st_ISGW2):
    _name = "FFBuToD2st@ISGW2"
    def __init__(self):
        super().__init__("B+", "D2*0")


@FFregistry.register
class LLSW(FFBToDstst.FFBToD2st_LLSW):
    _name = "FFBuToD2st@LLSW"
    def __init__(self):
        super().__init__("B+", "D2*0")


@FFregistry.register
class BLR(FFBToDstst.FFBToD2st_BLR):
    _name = "FFBuToD2st@BLR"
    def __init__(self):
        super().__init__("B+", "D2*0")