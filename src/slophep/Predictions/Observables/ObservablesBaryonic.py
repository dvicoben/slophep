from slophep.Predictions.FormFactorsBaryonic import FormFactorOneHalfpToOneHalfp
from slophep.Predictions.Observables import LbToOneHalfpEllNuPrediction



class LbToLcEllNuPrediction(LbToOneHalfpEllNuPrediction):
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorOneHalfpToOneHalfp,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("Lambdab", "Lambdac", "bc", lep, nu, FF, ffargs, par, scale)
