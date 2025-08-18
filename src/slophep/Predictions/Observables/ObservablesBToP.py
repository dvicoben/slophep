from slophep.Predictions.FormFactorsBToP import FormFactorBToP
from slophep.Predictions.Observables import BToPEllNuPrediction

class BpToDEllNuPrediction(BToPEllNuPrediction):
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("B+", "D0", "bc", lep, nu, FF, ffargs, par, scale)


class BdToDEllNuPrediction(BToPEllNuPrediction):
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("B0", "D+", "bc", lep, nu, FF, ffargs, par, scale)