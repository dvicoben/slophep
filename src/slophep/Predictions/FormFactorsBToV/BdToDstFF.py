import slophep.Predictions.FormFactorsBToV as FFBToV

class BLPR(FFBToV.BLPR_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BLPR"


class CLN(FFBToV.CLN_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_CLN"


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


class BGL_Hammer(FFBToV.BGL_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BGL_Hammer"
        self._ffpar = {
            "a0" : 0.00038,
            "a1" : 0.026905,
            "a2" : 0.,
            "b0" : 0.00055,
            "b1" : -0.0020370,
            "b2" : 0.,
            "c1" : -0.000433,
            "c2" : 0.005353,
            "d0" : 0.007,
            "d1" : -0.036
        }
    
    def get_ff(self, q2: float):
        etaEWVcb = self.internalparams["etaEW"]*self.internalparams["Vcb"]
        ff = super().get_ff(q2)
        return {k : ff[k]/etaEWVcb for k in ff}


class BSZ(FFBToV.BSZ_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_BSZ"


class HPQCD(FFBToV.HPQCD_BToV):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D*+", par, scale, *ffargs)
        self._name = "BdToDst_HPQCD"