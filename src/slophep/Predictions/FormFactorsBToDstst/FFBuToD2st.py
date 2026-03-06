import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD2st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D2*0", par, scale, *ffargs)
        self._name = "BuToD2st_ISGW2"