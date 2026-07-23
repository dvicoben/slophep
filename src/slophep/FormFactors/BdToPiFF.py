from typing import Any
import slophep.FormFactors.FormFactorsBToP as FFBToP

class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBdToPi@BSZ"
    def __init__(self):
        super().__init__("B0", "pi+")

    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "f+_0" : 0.1959500850585235,
            "f+_1" : -0.7763565146997043,
            "f+_2" : -0.19660537219243024,
            "f0_1" : -0.01446803146697799,
            "f0_2" : -0.19529667411047308,
            "fT_0" : 0.17311654822482725,
            "fT_1" : -0.40903716524130396,
            "fT_2" : -0.2580355775883392,
            # internalparams
            "m0": 5.540,
            "mp": 5.325
        }
        return ffpar


class BCL(FFBToP.FFBToP_BCL):
    _name = "FFBdToPi@BCL"
    def __init__(self):
        super().__init__("B0", "pi+")