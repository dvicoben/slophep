from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD0st
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1st
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD2st
from slophep.Predictions.Observables.BToD1EllNuObs import BToD1EllNuPrediction
from slophep.Predictions.Observables.BToD0stEllNuObs import BToD0stEllNuPrediction
from slophep.Predictions.Observables.BToD1stEllNuObs import BToD1stEllNuPrediction
from slophep.Predictions.Observables.BToD2stEllNuObs import BToD2stEllNuPrediction

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


class BuToD0stEllNuPrediction(BToD0stEllNuPrediction):
   def __init__(self, 
                lep: str, 
                nu: str,
                FF: FormFactorBToD0st,
                ffargs: list = [],
                par: dict = None,
                scale: float = 4.8,
                ):
        super().__init__("B+", "D0*0", "bc", lep, nu, FF, ffargs, par, scale)


class BuToD1stEllNuPrediction(BToD1stEllNuPrediction):
   def __init__(self, 
                lep: str, 
                nu: str,
                FF: FormFactorBToD1st,
                ffargs: list = [],
                par: dict = None,
                scale: float = 4.8,
                ):
        super().__init__("B+", "D1*0", "bc", lep, nu, FF, ffargs, par, scale)


class BuToD2stEllNuPrediction(BToD2stEllNuPrediction):
   def __init__(self, 
                lep: str, 
                nu: str,
                FF: FormFactorBToD2st,
                ffargs: list = [],
                par: dict = None,
                scale: float = 4.8,
                ):
        super().__init__("B+", "D2*0", "bc", lep, nu, FF, ffargs, par, scale)