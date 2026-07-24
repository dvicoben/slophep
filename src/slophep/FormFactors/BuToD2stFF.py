import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD2st_ISGW2):
    _name = "BuToD2st@ISGW2"
    def __init__(self):
        super().__init__("B+", "D2*0")


class LLSW(FFBToDstst.FFBToD2st_LLSW):
    _name = "BuToD2st@LLSW"
    def __init__(self):
        super().__init__("B+", "D2*0")


class BLR(FFBToDstst.FFBToD2st_BLR):
    _name = "BuToD2st@BLR"
    def __init__(self):
        super().__init__("B+", "D2*0")