import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)

class BLPR(FFBToP.BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)

class BGL(FFBToP.BGL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)

class CLN(FFBToP.CLN_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)