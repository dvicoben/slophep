"""
B0->D2*+ Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD2st_ISGW2):
    _name = "BdToD2st@ISGW2"
    def __init__(self):
        super().__init__("B0", "D2*+")


class LLSW(FFBToDstst.FFBToD2st_LLSW):
    _name = "BdToD2st@LLSW"
    def __init__(self):
        super().__init__("B0", "D2*+")


class BLR(FFBToDstst.FFBToD2st_BLR):
    _name = "BdToD2st@BLR"
    def __init__(self):
        super().__init__("B0", "D2*+")