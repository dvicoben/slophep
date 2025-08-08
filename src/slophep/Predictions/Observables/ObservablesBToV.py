from slophep.Predictions.FormFactorsBToV import FormFactorBToV
from slophep.Predictions.Observables import BToVEllNuPrediction

class BToDstEllNuPrediction(BToVEllNuPrediction):
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("B0", "D*+", "bc", lep, nu, FF, ffargs, par, scale)
        


class BsToDsstEllNuPrediction(BToVEllNuPrediction):
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("Bs", "Ds*", "bc", lep, nu, FF, ffargs, par, scale)



