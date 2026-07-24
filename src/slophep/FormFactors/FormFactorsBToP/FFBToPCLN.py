import numpy as np
from slophep.FormFactors.FormFactorsBToP.FFBToPBase import FormFactorBToP
from slophep.Tools.SamplingTools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)

class FFBToP_CLN(FormFactorBToP):
    _name = "FFBToP@CLN"
    def __init__(self, B: str, P: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, P)

    def define_userparams(self):
        ffpar = {
            "RhoSq"  : 1.186,
            "G1"     : 1.082,
            "Delta"  : 1.0  ,
            "a"      : 1.0  ,  # zero recoil expansion
            "mcOnmb" : 0.29 ,
            "lamBar" : 0.48 ,
            "ash"    : 0.26/np.pi
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Calculates CLN FFs.
        Implementation lifted from Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDCLN.cc, 

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
        rC = Mc/Mb
        sqrt2 = np.sqrt(2.)

        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1)
        a = self.get_userparam("a")
        zCon = (np.sqrt(w+1) - sqrt2*a)/(np.sqrt(w+1) + sqrt2*a)
        
        RhoSq = self.get_userparam("RhoSq")
        G1par = self.get_userparam("G1")
        Delta = self.get_userparam("Delta")

        V1wOnV1 = 1. - 8.* RhoSq * zCon + (51. * RhoSq - 10.) * zCon * zCon - (252. * RhoSq - 84.) * zCon * zCon *zCon
        fp=((1. + rC)/(2.*np.sqrt(rC))*V1wOnV1)*G1par
        S1wOnV1w = 1 + Delta*(-0.019 + 0.041*(w - 1.) - 0.015*np.power(w - 1.,2.))
        f0 = (np.sqrt(rC)/(1. + rC)*(1+w)*S1wOnV1w*V1wOnV1)*G1par

        ff = {
            "f+" : fp,
            "f0" : f0,
            "fT" : 0.0
        }
        return ff