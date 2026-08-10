"""
Lb->Lc Form-factors
"""
import slophep.FormFactors.FormFactorsBaryonic as FFBar
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class DKMR(FFBar.FFBaryonic_DKMR):
    _name = "FFLbToLc@DKMR"
    def __init__(self):
        super().__init__("Lambdab", "Lambdac")