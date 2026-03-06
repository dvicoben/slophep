import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD0st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0*0", par, scale, *ffargs)
        self._name = "BuToD0st_ISGW2"


class BLR(FFBToDstst.BLR_BToD0st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0*0", par, scale, *ffargs)
        self._name = "BuToD0st_BLR"