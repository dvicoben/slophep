"""
Bc->J/psi Form-factors
"""
from typing import Any
import numpy as np

# import slophep.FormFactors.FormFactorsBToV as FFBToV
from slophep.FormFactors.FormFactorsBToV import FormFactorBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from slophep.Core.user_registry import FFregistry

import logging
logger = logging.getLogger(__name__)

@FFregistry.register
class HPQCD2020(FormFactorBToV):
    _name = "FFBcToJpsi@HPQCD2020"
    def __init__(self):
        super().__init__("Bc", "J/psi")
        logger.info(f"{self.name} tensor FFs are zero.")

    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "a^0_A0" : 0.100612441741  ,
            "a^0_A1" : 0.0553196797329 ,
            "a^0_A2" : 0.0511266359181 ,
            "a^0_V"  : 0.105708947439  ,
            "a^1_A0" : -0.731340687341 ,
            "a^1_A1" : -0.265802751469 ,
            "a^1_A2" : -0.218722041876 ,
            "a^1_V"  : -0.74602264783  ,
            "a^2_A0" : 0.297215089789  ,
            "a^2_A1" : 0.310721659816  ,
            "a^2_A2" : -0.360026759263 ,
            "a^2_V"  : 0.102706760578  ,
            "a^3_A0" : -0.0220491268288,
            "a^3_A1" : 0.10616176997   ,
            "a^3_A2" : -0.0524140146237,
            "a^3_V"  : 0.00646710273132,
            "ResA0"  : np.array([6.2749, 6.872, 7.25])  ,
            "ResA1"  : np.array([6.745,6.75,7.15,7.15]) ,
            "ResA2"  : np.array([6.745,6.75,7.15,7.15]) ,
            "ResV"   : np.array([6.335,6.926,7.02,7.28]),
        }
        return ffpar

    def get_coef_arr(self, ff: str, nmax: int) -> np.ndarray:
        return np.array([self.get_userparam(f"a^{i}_{ff}" for i in range(nmax))])

    def pole(self, q2: float, tp: float, mres: np.ndarray) -> float:
        poles = (np.sqrt(tp - q2) - np.sqrt(tp - mres**2)) / (np.sqrt(tp - q2) + np.sqrt(tp - mres**2))
        return np.prod(poles) 

    def t0(self, tp: float, tm: float) -> float:
        return tm

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Calculate form-factors following BGL-like continuum fit in https://arxiv.org/pdf/2007.06957, Eq. (38). 
        Computation follows supplementary material `load_fit.py`.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict[str, float]
            Form-factor values
        """
        nmax = 4
        mB = self.get_param(f"m_{self.B}")
        mV = self.get_param(f"m_{self.V}")
        tp = (mB+mV)**2
        tm = (mB-mV)**2
        t0 = self.t0(tp, tm)
        z = (np.sqrt(tp - q2) - np.sqrt(tp - t0)) / (np.sqrt(tp - q2) + np.sqrt(tp - t0))
        zvec = np.array([z**k for k in range(nmax)])

        A0 = np.dot(zvec, self.get_coef_arr("A0", nmax))/self.pole(q2, tp, self.get_userparam("ResA0"))
        A1 = np.dot(zvec, self.get_coef_arr("A1", nmax))/self.pole(q2, tp, self.get_userparam("ResA1"))
        A2 = np.dot(zvec, self.get_coef_arr("A2", nmax))/self.pole(q2, tp, self.get_userparam("ResA2"))
        V  = np.dot(zvec, self.get_coef_arr("V" , nmax))/self.pole(q2, tp, self.get_userparam("ResV" ))
        A12 = ((A1 * (mB + mV)**2 * (mB**2 - mV**2 - q2)- A2 * (mB**4 + (mV**2 - q2)**2
                - 2 * mB**2 * (mV**2 + q2))) / (16. * mB * mV**2 * (mB + mV)))

        ff = {
            "A0"  : A0,
            "A1"  : A1,
            "A2"  : A2,
            "A12" : A12,
            "V"   : V,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "T3"  : 0.0,
            "T23" : 0.0
        }
        return ff