import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BSZ"


class BLPR(FFBToP.BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BLPR"


class BLPR_Hammer(FFBToP.BLPR_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BLPR_Hammer"


class BLPRXP(FFBToP.BLPRXP_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BLPRXP"


class BGL(FFBToP.BGL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BGL"


class BGL_Hammer(FFBToP.BGL_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BGL_Hammer"


class CLN(FFBToP.CLN_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_CLN"


class CLN_Hammer(FFBToP.CLN_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_CLN_Hammer"