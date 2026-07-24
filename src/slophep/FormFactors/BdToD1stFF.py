import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD1st_ISGW2):
    _name = "BdToD1st@ISGW2"
    def __init__(self):
        super().__init__("B0", "D1*+")


class LLSW(FFBToDstst.FFBToD1st_LLSW):
    _name = "BdToD1st@LLSW"
    def __init__(self):
        super().__init__("B0", "D1*+")


class BLR(FFBToDstst.FFBToD1st_BLR):
    _name = "BdToD1st@BLR"
    def __init__(self):
        super().__init__("B0", "D1*+")