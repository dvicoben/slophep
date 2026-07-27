"""
B+->D1 Form-factors
"""
import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD1_ISGW2):
    _name = "BuToD1@ISGW2"
    def __init__(self):
        super().__init__("B+", "D10")


class LLSW(FFBToDstst.FFBToD1_LLSW):
    _name = "BuToD1@LLSW"
    def __init__(self):
        super().__init__("B+", "D10")


class BLR(FFBToDstst.FFBToD1_BLR):
    _name = "BuToD1@BLR"
    def __init__(self):
        super().__init__("B+", "D10")