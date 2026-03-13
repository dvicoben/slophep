import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class BLR(FFBToDstst.BLR_BToD1):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D1+", par, scale, *ffargs)
        self._name = "BdToD1_BLR"


class ISGW2(FFBToDstst.ISGW2_BToD1):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D1+", par, scale, *ffargs)
        self._name = "BdToD1_ISGW2"