import slophep.Predictions.FormFactorsBToDstst as FFBToDstst

class ISGW2(FFBToDstst.ISGW2_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D1*+", par, scale, *ffargs)
        self._name = "BdToD1st_ISGW2"


class LLSW(FFBToDstst.LLSW_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D1*+", par, scale, *ffargs)
        self._name = "BdToD1st_LLSW"


class BLR(FFBToDstst.BLR_BToD1st):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D1*+", par, scale, *ffargs)
        self._name = "BdToD1st_BLR"