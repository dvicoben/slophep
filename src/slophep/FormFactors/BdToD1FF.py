"""
B0->D1 Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD1_ISGW2):
    _name = "BdToD1@ISGW2"
    def __init__(self):
        super().__init__("B0", "D1+")


@FFregistry.register
class LLSW(FFBToDstst.FFBToD1_LLSW):
    _name = "BdToD1@LLSW"
    def __init__(self):
        super().__init__("B0", "D1+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD1_BLR):
    _name = "BdToD1@BLR"
    def __init__(self):
        super().__init__("B0", "D1+")