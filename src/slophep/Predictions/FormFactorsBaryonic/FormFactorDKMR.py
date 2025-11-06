import numpy as np
from math import sqrt

from slophep.Predictions.FormFactorsBaryonic.FormFactorOneHalfpToOneHalfp import FormFactorOneHalfpToOneHalfp


class DKMR_OneHalfpToOneHalfp(FormFactorOneHalfpToOneHalfp):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, C, par, scale)

        self._name = "OneHalfpToOneHalfp_DKMR"
        self._ffpar = {
            "aV_t0"  : 0.74392566003686,
            "aV_t1"  : -4.6476511382467,
            "aV_t2"  : 0.0,
            "aA_t0"  : 0.73960499809758,
            "aA_t1"  : -4.3664554520299,
            "aA_t2"  : 0.0,
            "aV_l0"  : 0.81458549275845,
            "aV_l1"  : -4.8987680713973,
            "aV_l2"  : 0.0,
            "aV_p0"  : 1.077996185173,
            "aV_p1"  : -6.4170836206134,
            "aV_p2"  : 0.0,
            "aA_l0"  : 0.6846557018386,
            "aA_l1"  : -4.4311552783827,
            "aA_l2"  : 0.0,
            "aA_p1"  : -4.4633624514401,
            "aA_p2"  : 0.0,
            "aT_l0"  : 0.97518189850197,
            "aT_l1"  : -5.4999842229709,
            "aT_l2"  : 0.0,
            "aT_p0"  : 0.70539533108116,
            "aT_p1"  : -4.3577976693726,
            "aT_p2"  : 0.0,
            "aT5_l0" : 0.67275080357035,
            "aT5_l1" : -4.4321757873186,
            "aT5_l2" : 0.0,
            "aT5_p1" : -4.4927847697461,
            "aT5_p2" : 0.0
        }
        self._params = [
            "aV_t0", "aV_t1", "aV_t2", "aA_t0", "aA_t1", "aA_t2", 
            "aV_l0", "aV_l1", "aV_l2", "aV_p0", "aV_p1", "aV_p2",
            "aA_l0", "aA_l1", "aA_l2", "aA_p1", "aA_p2", 
            "aT_l0", "aT_l1", "aT_l2", "aT_p0", "aT_p1", "aT_p2",
            "aT5_l0", "aT5_l1", "aT5_l2", "aT5_p1", "aT5_p2"
        ]
        # Resonance masses from https://arxiv.org/pdf/1702.02243 (Table 1) and https://arxiv.org/pdf/1503.01421 (Table 7)
        internalparams = {
            "m1p" : 6.768,
            "m1m" : 6.332,
            "m0p" : 6.725,
            "m0m" : 6.276
        }
        self._internalparams.update(internalparams)

    def z(self, q2: float, mres: float):
        # eqns 3.3-3.5 in https://arxiv.org/pdf/1702.02243
        mB = self.internalparams["Mb"]
        mM = self.internalparams["Mc"]
        t0 = (mB-mM)**2
        tp = (mres)**2
        sq2 = sqrt(tp-q2)
        st0 = sqrt(tp-t0)
        return (sq2-st0)/(sq2+st0)
    
    def get_ff(self, q2: float) -> dict[str, float]:
        
        m0p = self.internalparams["m0p"]
        z0p = self.z(q2, m0p)
        m0m = self.internalparams["m0m"]
        z0m = self.z(q2, m0m)
        m1p = self.internalparams["m1p"]
        z1p = self.z(q2, m1p)
        m1m = self.internalparams["m1m"]
        z1m = self.z(q2, m1m)

        fV_t = 1.0/(1. - (q2/(m0p**2)))*(
            self.ffpar["aV_t0"]
            + self.ffpar["aV_t1"]*z0p
            + self.ffpar["aV_t2"]*(z0p**2)
        )
        fV_l = 1.0/(1. - (q2/(m1m**2)))*(
            self.ffpar["aV_l0"]
            + self.ffpar["aV_l1"]*z1m
            + self.ffpar["aV_l2"]*(z1m**2)
        )
        fV_p = 1.0/(1. - (q2/(m1m**2)))*(
            self.ffpar["aV_p0"]
            + self.ffpar["aV_p1"]*z1m
            + self.ffpar["aV_p2"]*(z1m**2)
        )
        fA_t = 1.0/(1. - (q2/(m0m**2)))*(
            self.ffpar["aA_t0"]
            + self.ffpar["aA_t1"]*z0m
            + self.ffpar["aA_t2"]*(z0m**2)
        )
        fA_l = 1.0/(1. - (q2/(m1p**2)))*(
            self.ffpar["aA_l0"]
            + self.ffpar["aA_l1"]*z1p
            + self.ffpar["aA_l2"]*(z1p**2)
        )
        fA_p = 1.0/(1. - (q2/(m1p**2)))*(
            self.ffpar["aA_l0"]
            + self.ffpar["aA_p1"]*z1p
            + self.ffpar["aA_p2"]*(z1p**2)
        )
        fT_l = 1.0/(1. - (q2/(m1m**2)))*(
            self.ffpar["aT_l0"]
            + self.ffpar["aT_l1"]*z1m
            + self.ffpar["aT_l2"]*(z1m**2)
        )
        fT_p = 1.0/(1. - (q2/(m1m**2)))*(
            self.ffpar["aT_p0"]
            + self.ffpar["aT_p1"]*z1m
            + self.ffpar["aT_p2"]*(z1m**2)
        )
        fT5_l = 1.0/(1. - (q2/(m1p**2)))*(
            self.ffpar["aT5_l0"]
            + self.ffpar["aT5_l1"]*z1p
            + self.ffpar["aT5_l2"]*(z1p**2)
        )
        fT5_p = 1.0/(1. - (q2/(m1p**2)))*(
            self.ffpar["aT5_l0"]
            + self.ffpar["aT5_p1"]*z1p
            + self.ffpar["aT5_p2"]*(z1p**2)
        )

        ff = {
            "Vt"  : fV_t,
            "V0"  : fV_l,
            "Vp"  : fV_p,
            "At"  : fA_t,
            "A0"  : fA_l,
            "Ap"  : fA_p,
            "T0"  : fT_l,
            "T50" : fT5_l,
            "Tp"  : fT_p,
            "T5p" : fT5_p
        }

        return ff