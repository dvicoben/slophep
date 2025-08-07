from slophep.Predictions.FormFactorsBToP.BToPFFBSZ import BSZ_BToP

class BSZ(BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)