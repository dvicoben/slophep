"""
B+->D0*0 Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD0st_ISGW2):
    _name = "BuToD0st@ISGW2"
    def __init__(self):
        super().__init__("B+", "D0*0")


@FFregistry.register
class LLSW(FFBToDstst.FFBToD0st_LLSW):
    _name = "BuToD0st@LLSW"
    def __init__(self):
        super().__init__("B+", "D0*0")


@FFregistry.register
class BLR(FFBToDstst.FFBToD0st_BLR):
    _name = "BuToD0st@BLR"
    def __init__(self):
        super().__init__("B+", "D0*0")