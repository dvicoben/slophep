from slophep.Observables.ObservableImpl import BToPEllNuPrediction, MixinBToPAngular
from slophep.FormFactors.FormFactorsBToP import FormFactorBToP


class BuToDEllNuPrediction(BToPEllNuPrediction, MixinBToPAngular):
    _name = "ObsBuToDEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ):
        super().__init__("B+", "D0", "bc", lep, nu, FF)



class BdToDEllNuPrediction(BToPEllNuPrediction, MixinBToPAngular):
    _name = "ObsBdToDEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ):
        super().__init__("B0", "D+", "bc", lep, nu, FF)



class BdToPiEllNuPrediction(BToPEllNuPrediction, MixinBToPAngular):
    _name = "ObsBdToPiEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ):
        super().__init__("B0", "pi+", "bu", lep, nu, FF)



class BsToKEllNuPrediction(BToPEllNuPrediction):
    _name = "ObsBsToKEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToP,
                 ):
        super().__init__("Bs", "K+", "bu", lep, nu, FF)
