"""
Bs->K Form-factors
"""
from typing import Any
import slophep.FormFactors.FormFactorsBToP as FFBToP
from slophep.Tools.errfluct_tools import fluctsettings, FluctType


class BCL(FFBToP.FFBToP_BCLGeneric):
    _name = "FFBsToK@BCL"
    def __init__(self):
        super().__init__("Bs", "K+", 4, 3)

    def define_userparams(self):
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f+_0" : 0.374 , 
            "a_f+_1" : -0.672, 
            "a_f+_2" : 0.07  , 
            "a_f+_3" : 1.34  ,
            "a_f0_0" : 0.2203, 
            "a_f0_1" : 0.089 , 
            "a_f0_2" : 0.24  ,
            #internalparams
            "m1m"    : 5.325,
            "m0p"    : 5.68 ,
            "tp"     : (5.280 + 0.13957)**2,
            "q2cons" : True, # impose f_+(q^2=0) = f_0(q^2=0)
        })
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        ff =  super().get_ff(q2)
        m0p = self.get_userparam("m0p")
        P0p = (1. - q2/(m0p*m0p))
        ffs = {
            "f+" : ff["f+"],
            "f0" : ff["f0"]*(1./P0p),
            "fT" : ff["fT"]
        }
        return ffs