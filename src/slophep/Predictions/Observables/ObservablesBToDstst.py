from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1
from slophep.Predictions.Observables.BToD1EllNuObs import BToD1EllNuPrediction

class BuToD1EllNuPrediction(BToD1EllNuPrediction):
   def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__("B+", "D10", "bc", lep, nu, FF, ffargs, par, scale)