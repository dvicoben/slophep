import slophep.Predictions.FormFactorsBToV as FFBToV


class ISGW2(FFBToV.ISGW2_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_ISGW2"


class BLPR(FFBToV.BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPR"


class BLPR_Hammer(FFBToV.BLPR_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPR_Hammer"


class BLPRXP(FFBToV.BLPRXP_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPRXP"


class CLN(FFBToV.CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN"


class CLN_Hammer(FFBToV.CLN_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN_Hammer"


class CLN2(FFBToV.CLN2_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN2"
        self.internalparams.update({
            "qiqj" : "bc"
        })


class BGL(FFBToV.BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BGL"


class BGL_Hammer(FFBToV.BGL_BToV_Hammer):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BGL_Hammer"


class BSZ(FFBToV.BSZ_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BSZ"


class HPQCD(FFBToV.HPQCD_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_HPQCD"