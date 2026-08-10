from slophep.Observables.ObservableImpl import BToVEllNuPrediction, MixinBToVAngular
from slophep.FormFactors.FormFactorsBToV import FormFactorBToV


class BdToDstEllNuPrediction(BToVEllNuPrediction, MixinBToVAngular):
    _name = "ObsBdToDstEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ):
        super().__init__("B0", "D*+", "bc", lep, nu, FF)
        


class BsToDsstEllNuPrediction(BToVEllNuPrediction, MixinBToVAngular):
    _name = "ObsBsToDsstEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ):
        super().__init__("Bs", "Ds*", "bc", lep, nu, FF)



class BcToJpsiEllNuPrediction(BToVEllNuPrediction):
    _name = "ObsBcToJpsiEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToV,
                 ):
        super().__init__("Bc", "J/psi", "bc", lep, nu, FF)
