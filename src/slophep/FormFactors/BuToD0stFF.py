import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD0st_ISGW2):
    _name = "BuToD0st@ISGW2"
    def __init__(self):
        super().__init__("B+", "D0*0")


class LLSW(FFBToDstst.FFBToD0st_LLSW):
    _name = "BuToD0st@LLSW"
    def __init__(self):
        super().__init__("B+", "D0*0")


class BLR(FFBToDstst.FFBToD0st_BLR):
    _name = "BuToD0st@BLR"
    def __init__(self):
        super().__init__("B+", "D0*0")