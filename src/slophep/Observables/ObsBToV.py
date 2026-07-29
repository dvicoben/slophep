from slophep.Observables.ObservableImpl import BToVEllNuPrediction
from slophep.FormFactors.FormFactorsBToV import FormFactorBToV


class BdToDstEllNuPrediction(BToVEllNuPrediction):
    _name = "ObsBdToDstEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ):
        super().__init__("B0", "D*+", "bc", lep, nu, FF)
        


class BsToDsstEllNuPrediction(BToVEllNuPrediction):
    _name = "ObsBsToDsstEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ):
        super().__init__("Bs", "Ds*", "bc", lep, nu, FF)



