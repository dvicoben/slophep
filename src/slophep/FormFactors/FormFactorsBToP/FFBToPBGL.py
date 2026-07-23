import numpy as np
from slophep.FormFactors.FormFactorsBToP.FFBToPBase import FormFactorBToP
from slophep.Tools.SamplingTools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


class FFBToP_BGLGeneric(FormFactorBToP):
    _name = "FFBToP@BGLGen"
    def __init__(self, B: str, P: str,
                 N_fp : int = 3,
                 N_f0 : int = 3):
        logger.warning(f"{self.name} tensor FFs are zero.")
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
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        }
        ffpar.update({f"a_f+_{iord}" : 0.0 for iord in range(self.n["f+"])})
        ffpar.update({f"a_f0_{iord}" : 0.0 for iord in range(self.n["f0"])})
        # ffpar.update({f"a_fT_{iord}" : 0.0 for iord in range(self.n["fT"])})
        return ffpar

    def get_coef_arr(self, ffstr: str) -> list[float]:
        return np.array([self.get_userparam(f"a_{ffstr}_{iord}") for iord in range(self.n[ffstr])])

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = np.sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr)

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
        Mc = self.get_param(f"m_{self.P}")
        Mb2 = Mb*Mb
        rC = Mc/Mb
        rC2 = rC*rC
        sqrC = np.sqrt(rC)

        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))
        zpow = np.array([z**ik for ik in range(int(self.internalparams["nmax"]))])

        ap = self.get_coef_arr("f+")
        a0 = self.get_coef_arr("f0")
        
        chiT = self.get_userparam("ChiT")
        chiL = self.get_userparam("ChiL")
        nc = self.get_userparam("nc")

        Pp = self.blaschke(self.get_userparam("BcStatesp"), z, Mb, Mc)
        P0 = self.blaschke(self.get_userparam("BcStates0"), z, Mb, Mc)

        phip = 32.*np.sqrt(nc/(6.*np.pi*chiT*Mb2))*rC2*pow(1+z,2)*pow(1-z,0.5)/pow((1+rC)*(1-z)+2*sqrC*(1+z), 5)
        phi0 = (1-rC2)*8.*np.sqrt(nc/(8.*np.pi*chiL))*rC*(1+z)*pow(1-z,1.5)/pow((1+rC)*(1-z)+2*sqrC*(1+z), 4)

        fp = np.dot(ap, zpow)/(Pp*phip)
        f0 = np.dot(a0, zpow)/(P0*phi0)

        ff = {
            "f+" : fp,
            "f0" : f0,
            "fT" : 0.0
        }
        return ff




class FFBToP_BGL(FFBToP_BGLGeneric):
    _name = "FFBToP@BGL"
    def __init__(self, B: str, P: str):
        logger.warning(f"{self.name} tensor FFs are zero.")
        super().__init__(B, P, 4, 4)

    def define_userparams(self):
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f+_0" : 0.01565,
            "a_f+_1" : -0.0353,
            "a_f+_2" : -0.043,
            "a_f+_3" : 0.194,
            "a_f0_0" : 0.07932,
            "a_f0_1" : -0.214,
            "a_f0_2" : 0.17,
            "a_f0_3" : -0.958,
            #internalparams
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        })
        return ffpar
