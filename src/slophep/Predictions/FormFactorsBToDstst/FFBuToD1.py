import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class BLR(FFBToDstst.BLR_BToD1):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D10", par, scale, *ffargs)
        self._name = "BuToD1_BLR"


class ISGW2(FFBToDstst.ISGW2_BToD1):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D10", par, scale, *ffargs)
        self._name = "BuToD1_ISGW2"