"""
B0->D0*+ Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD0st_ISGW2):
    _name = "FFBdToD0st@ISGW2"
    def __init__(self):
        super().__init__("B0", "D0*+")


@FFregistry.register
class LLSW(FFBToDstst.FFBToD0st_LLSW):
    _name = "FFBdToD0st@LLSW"
    def __init__(self):
        super().__init__("B0", "D0*+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD0st_BLR):
    _name = "FFBdToD0st@BLR"
    def __init__(self):
        super().__init__("B0", "D0*+")