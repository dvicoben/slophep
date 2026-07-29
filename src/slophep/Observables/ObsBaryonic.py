from slophep.Observables.ObservableImpl import LbToOneHalfpEllNuPrediction
from slophep.FormFactors.FormFactorsBaryonic import FormFactorOneHalfpToOneHalfp


class LbToLcEllNuPrediction(LbToOneHalfpEllNuPrediction):
    _name = "ObsLbToLcEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorOneHalfpToOneHalfp,
                 ):
        super().__init__("Lambdab", "Lambdac", "bc", lep, nu, FF)

