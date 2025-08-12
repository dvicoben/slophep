from slophep.Predictions.FormFactorsBToP import BSZ_BToP
from slophep.Predictions.FormFactorsBToP import BLPR_BToP
from slophep.Predictions.FormFactorsBToP import BGL_BToP
from slophep.Predictions.FormFactorsBToP import CLN_BToP

class BSZ(BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)

class BLPR(BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)

class BGL(BGL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)

class CLN(CLN_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)