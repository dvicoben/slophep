import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D2*+", par, scale, *ffargs)
        self._name = "BdToD2st_ISGW2"


class BLR(FFBToDstst.BLR_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D2*+", par, scale, *ffargs)
        self._name = "BdToD2st_BLR"