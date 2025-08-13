import slophep.Predictions.FormFactorsBToV as FFBToV

class BLPR(FFBToV.BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class CLN(FFBToV.CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class CLN2(FFBToV.CLN2_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class BGL(FFBToV.BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class BSZ(FFBToV.BSZ_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class HPQCD(FFBToV.HPQCD_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)