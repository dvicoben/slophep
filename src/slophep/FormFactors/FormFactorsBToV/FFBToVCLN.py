import numpy as np
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.SamplingTools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)

class FFBToV_CLN(FormFactorBToV):
    _name = "FFBToV@CLN"
    
    def __init__(self, B: str, V: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, V)

    def define_userparams(self):
        ffpar = {
            "RhoSq" : 1.207,
            "h_A1"  : 0.908,
            "R1"    : 1.401,
            "R2"    : 0.854,
            "R0"    : 1.15
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Central value of $B\to V$ form factors in the lattice convention
        CLN parametrization. See eqs. (B4)-(B7) of arXiv:1203.2654.

        Directly lifted from flavio, see https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py

        For SM, this behaviour should be the same as HQET2 in EvtGen and CLN in HAMMER (though haven't verified throughly)

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict
            FF dictionary
        """

        mB = self.get_param(f"m_{self.B}")
        mV = self.get_param(f"m_{self.V}")
        w = (mB**2 + mV**2 - q2) / (2*mB*mV)
        z = (np.sqrt(w+1)-np.sqrt(2))/(np.sqrt(w+1)+np.sqrt(2))
        RV = 2*np.sqrt(mB*mV)/(mB+mV)

        hA1_1 = self.get_userparam('h_A1')
        R1_1  = self.get_userparam('R1')
        R2_1  = self.get_userparam('R2')
        R0_1  = self.get_userparam('R0')
        rho2  = self.get_userparam('RhoSq')

        hA1 = hA1_1 * (1 - 8*rho2*z + (53*rho2-15)*z**2 - (231*rho2-91)*z**3)
        R1 = R1_1 - 0.12*(w-1) + 0.05*(w-1)**2
        R2 = R2_1 + 0.11*(w-1) - 0.06*(w-1)**2
        # R0 from Fajfer et al, Phys. Rev. D 85 094025 (2012)
        R0 = R0_1 - 0.11*(w-1) + 0.01*(w-1)**2

        ff = {}
        ff['A1'] = hA1 * RV * (w+1)/2.
        ff['A0'] = R0/RV * hA1
        # A2 is not used in the lattice convention
        A2 = R2/RV * hA1
        # conversion from A_1, A_2 to A_12
        ff['A12'] = ((ff['A1']*(mB + mV)**2* (mB**2 - mV**2 - q2)
                - A2*(mB**4 + (mV**2 - q2)**2 - 2*mB**2*(mV**2 + q2)))
                / (16.*mB*mV**2*(mB + mV)))
        ff['V'] = R1/RV * hA1
        # SM only non tensor FFs
        ff["T1"] = 0.0
        ff["T2"] = 0.0
        ff["T23"] = 0.0
        return ff