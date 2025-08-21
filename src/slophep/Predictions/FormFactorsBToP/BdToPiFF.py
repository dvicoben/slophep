import slophep.Predictions.FormFactorsBToP as FFBToP

class BSZ(FFBToP.BSZ_BToP):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__("B0", "pi+", par, scale, *ffargs)
        
        self._ffpar = {
            "f+_0" : 0.1959500850585235,
            "f+_1" : -0.7763565146997043,
            "f+_2" : -0.19660537219243024,
            "f0_1" : -0.01446803146697799,
            "f0_2" : -0.19529667411047308,
            "fT_0" : 0.17311654822482725,
            "fT_1" : -0.40903716524130396,
            "fT_2" : -0.2580355775883392
        }
        internalparams = {
            'm0': 5.540,
            'mp': 5.325,
        }
        self._internalparams.update(internalparams)