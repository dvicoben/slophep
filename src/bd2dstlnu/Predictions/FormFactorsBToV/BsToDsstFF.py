from bd2dstlnu.Predictions.FormFactorsBToV.BToVFFBLPR import BLPR_BToV
from bd2dstlnu.Predictions.FormFactorsBToV.BToVFFCLN import CLN_BToV
from bd2dstlnu.Predictions.FormFactorsBToV.BToVFFBGL import BGL_BToV
from bd2dstlnu.Predictions.FormFactorsBToV.BToVFFHPQCD import HPQCD_BToV



class BLPR(BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds*", par, scale, *ffargs)

class CLN(CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds*", par, scale, *ffargs)

class BGL(BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "Ds*", par, scale, *ffargs)