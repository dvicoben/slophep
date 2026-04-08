import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("Bs", "K+", par, scale, *ffargs)
        # Parameters from B->P EOS https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-p-form-factor-parametrizations
        # and https://arxiv.org/pdf/1811.00983
        self._ffpar = {
            "f+_0" : 0.32908992854636815 ,
            "f+_1" : -0.8669465867427438 ,
            "f+_2" : 0.006095669673341246,
            "f0_1" : 0.19511719957037596 ,
            "f0_2" : -0.4461264576740431 ,
            "fT_0" : 0.2993831492500477  ,
            "fT_1" : -0.7735456824025474 ,
            "fT_2" : 0.009554375102511888
        }
        internalparams = {
            'm0': 5.630,
            'mp': 5.412,
        }
        self._internalparams.update(internalparams)