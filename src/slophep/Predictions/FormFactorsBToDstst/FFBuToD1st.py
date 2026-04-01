import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D1*0", par, scale, *ffargs)
        self._name = "BuToD1st_ISGW2"


class LLSW(FFBToDstst.LLSW_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D1*0", par, scale, *ffargs)
        self._name = "BuToD1st_LLSW"


class BLR(FFBToDstst.BLR_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D1*0", par, scale, *ffargs)
        self._name = "BuToD1st_BLR"