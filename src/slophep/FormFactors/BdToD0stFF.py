import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD0st_ISGW2):
    _name = "BdToD0st@ISGW2"
    def __init__(self):
        super().__init__("B0", "D0*+")


class LLSW(FFBToDstst.FFBToD0st_LLSW):
    _name = "BdToD0st@LLSW"
    def __init__(self):
        super().__init__("B0", "D0*+")


class BLR(FFBToDstst.FFBToD0st_BLR):
    _name = "BdToD0st@BLR"
    def __init__(self):
        super().__init__("B0", "D0*+")