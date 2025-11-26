import slophep.Predictions.FormFactorsBaryonic as FFBar

class DKMR(FFBar.DKMR_OneHalfpToOneHalfp):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Lambdab", "Lambdac", par, scale, *ffargs)
        self._name = "LbToLc_DKMR"