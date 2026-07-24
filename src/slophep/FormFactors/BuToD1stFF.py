import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.FFBToD1st_ISGW2):
    _name = "BuToD1st@ISGW2"
    def __init__(self):
        super().__init__("B+", "D1*0")


class LLSW(FFBToDstst.FFBToD1st_LLSW):
    _name = "BuToD1st@LLSW"
    def __init__(self):
        super().__init__("B+", "D1*0")


class BLR(FFBToDstst.FFBToD1st_BLR):
    _name = "BuToD1st@BLR"
    def __init__(self):
        super().__init__("B+", "D1*0")