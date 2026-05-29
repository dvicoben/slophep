import slophep.Predictions.FormFactorsBToP as FFBToP


class BCL(FFBToP.BCL_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "K+", par, scale, *ffargs)
        self._ffpar = {
            "f+_0" : 0.374 , 
            "f+_1" : -0.672, 
            "f+_2" : 0.07  , 
            "f+_3" : 1.34  ,
            "f0_0" : 0.2203, 
            "f0_1" : 0.089 , 
            "f0_2" : 0.24
        }
        internalparams = {
            "m1m"    : 5.325,
            "m0p"    : 5.68 ,
            "tp"     : (5.280 + 0.13957)**2,
            "q2cons" : True, # impose f_+(q^2=0) = f_0(q^2=0)
        }
        self._internalparams.update(internalparams)
    
    def get_ff_mmeson(self, q2: float, mM: float, mB: float = None) -> dict[str, float]:
        # https://arxiv.org/pdf/2111.09849, https://arxiv.org/pdf/1901.02561
        ff =  super().get_ff_mmeson(q2, mM, mB)
        m0p = self.internalparams["m0p"]
        P0p = (1. - q2/(m0p*m0p))
        ffs = {
            "f+" : ff["f+"],
            "f0" : ff["f0"]*(1./P0p),
            "fT" : ff["fT"]
        }
        return ffs