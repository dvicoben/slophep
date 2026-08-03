import numpy as np
from slophep.FormFactors.FormFactorsBToP.FFBToPBase import FormFactorBToP
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


class FFBToP_BCLGeneric(FormFactorBToP):
    _name = "FFBToP@BCLGen"
    def __init__(self, B: str, P: str,
                 N_fp : int, N_f0 : int):
        logger.info(f"{self.name} tensor FFs are zero.")
        self._n = {
            "f+"   : N_fp,
            "f0"   : N_f0,
            "nmax" : max(N_fp, N_f0)
        }
        super().__init__(B, P)

    @property
    def n(self) -> dict[str, int]:
        return self._n

    def define_userparams(self):
        ffpar = {
            "m1m"    : 5.325,
            "q2cons" : False,   # impose f_+(q^2=0) = f_0(q^2=0)
            "tp"     : None
        }
        ffpar.update({f"a_f+_{iord}" : 0.0 for iord in range(self.n["f+"])})
        ffpar.update({f"a_f0_{iord}" : 0.0 for iord in range(self.n["f0"])})
        # ffpar.update({f"a_fT_{iord}" : 0.0 for iord in range(self.n["fT"])})
        return ffpar

    def get_coef_arr(self, ffstr: str) -> list[float]:
        return np.array([self.get_userparam(f"a_{ffstr}_{iord}") for iord in range(self.n[ffstr])])
    
    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict:
        """Calculates BGL FFs.
        Implementation lifted from Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBGL.cc, 

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        Mb = self.get_param(f"m_{self.B}")
        Mu = self.get_param(f"m_{self.P}")

        m1m = self.get_userparam("m1m")
        P1m = 1-(q2/(m1m*m1m))

        # w = (Mb**2 + Mu**2 - q2) / (2. * Mb * Mu)
        tp = self.get_userparam("tp") if self.get_userparam("tp") is not None else (Mb + Mu)*(Mb + Mu)
        tm = (Mb - Mu)*(Mb - Mu)
        t0 = tp*(1. - np.sqrt(1. - tm/tp))
        z = (np.sqrt(tp - q2) - np.sqrt(tp - t0))/(np.sqrt(tp - q2) + np.sqrt(tp - t0))
        z0 = (np.sqrt(tp) - np.sqrt(tp - t0))/(np.sqrt(tp) + np.sqrt(tp - t0))

        fp_vec = self.get_coef_arr("f+")
        f0_vec = self.get_coef_arr("f0")
        Nz = len(fp_vec)
        N0 = len(f0_vec)
        Nmax = max(Nz, N0)
        z_vec  = np.array([z**k for k in range(Nmax+1)])
        neg1_z = np.array([((-1.)**(n-Nz))*(n/Nz)*z_vec[Nz] for n in range(Nz)])

        f0 = np.sum(f0_vec*z_vec[:N0])
        fp = (1/P1m)*np.sum(
            fp_vec*(z_vec[:Nz] - neg1_z)
        )

        if self.get_userparam("q2cons"):
            z0_vec = np.array([z0**k for k in range(Nmax+1)])
            neg1_z0 = np.array([(-1.)**(n-Nz)*(n/Nz)*z0_vec[Nz] for n in range(Nz)])

            fpq2 = np.sum(fp_vec*(z0_vec[:Nz] - neg1_z0[:Nz]))
            f0q2 = np.sum(f0_vec*z0_vec[:N0])
            f0 += (fpq2 - f0q2)/z0_vec[N0]*z_vec[N0]

        ff = {
            "f+" : fp ,
            "f0" : f0 ,
            "fT" : 0.0
        }
        return ff



class FFBToP_BCL(FFBToP_BCLGeneric):
    _name = "FFBToP@BCL"
    def __init__(self, B: str, P: str):
        super().__init__(B, P, 4, 4)

    def define_userparams(self):
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f+_0" : 0.419 ,
            "a_f+_1" : -0.495,
            "a_f+_2" : -0.43 ,
            "a_f+_3" : 0.22  ,
            "a_f0_0" : 0.510 ,
            "a_f0_1" : -1.700,
            "a_f0_2" : 1.53  ,
            "a_f0_3" : 4.52  ,
            #internalparams
            "m1m"    : 5.325  ,
            "q2cons" : False   # impose f_+(q^2=0) = f_0(q^2=0)
        })
        return ffpar