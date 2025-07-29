from slophep.Predictions.FormFactorsBToV.BToVFFBLPR import BLPR_BToV
from slophep.Predictions.FormFactorsBToV.BToVFFCLN import CLN_BToV
from slophep.Predictions.FormFactorsBToV.BToVFFBGL import BGL_BToV
from slophep.Predictions.FormFactorsBToV.BToVFFBSZ import BSZ_BToV
from slophep.Predictions.FormFactorsBToV.BToVFFHPQCD import HPQCD_BToV



class BLPR(BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class CLN(CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class BGL(BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class BSZ(BSZ_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)

class HPQCD(HPQCD_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)