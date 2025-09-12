import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_BGL"


class BLPR(FFBToP.BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_BLPR"


class BLPR_Hammer(FFBToP.BLPR_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_BLPR_Hammer"


class BGL(FFBToP.BGL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_BGL"


class BGL_Hammer(FFBToP.BGL_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_BGL_Hammer"


class CLN(FFBToP.CLN_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_CLN"


class CLN_Hammer(FFBToP.CLN_BToP_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B+", "D0", par, scale, *ffargs)
        self._name = "BuToD_CLN_Hammer"