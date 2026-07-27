"""
Lb->Lc Form-factors
"""
import slophep.FormFactors.FormFactorsBaryonic as FFBar

class DKMR(FFBar.FFBaryonic_DKMR):
    _name = "FFLbToLc@DKMR"
    def __init__(self):
        super().__init__("Lambdab", "Lambdac")