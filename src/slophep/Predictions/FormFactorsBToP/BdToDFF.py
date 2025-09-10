import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BSZ"


class BLPR(FFBToP.BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BLPR"


class BLPR_Hammer(FFBToP.BLPR_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BLPR_Hammer"
        print(f"WARNING: {self.name} uses a different basis for fT than the one assumed by flavio to compute observables. "
              +"Take care if using this for anything other than FF comparisons!")
    
    def get_ff(self, q2: float) -> dict:
        ff = super().get_ff(q2)
        ff["fT"] = ff["fT"]/(self.internalparams["Mb"] + self.internalparams["Mc"])
        return ff


class BGL(FFBToP.BGL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_BGL"


class CLN(FFBToP.CLN_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "D+", par, scale, *ffargs)
        self._name = "BdToD_CLN"