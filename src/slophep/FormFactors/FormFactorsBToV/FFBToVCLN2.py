import numpy as np
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from flavio.physics.running import running

import logging
logger = logging.getLogger(__name__)

class FFBToV_CLN2(FormFactorBToV):
    """CLN form-factors, tensor FFs obtained using Eq. 11 in https://arxiv.org/pdf/1503.05534"""
    _name = "FFBToV@CLN2"
    
    def __init__(self, B: str, V: str):
        logger.info(f"{self.name} Tensor FFs are set assuming eq. 11 in https://arxiv.org/pdf/1503.05534 is unity, use with care for BSM.")
        super().__init__(B, V)

    def define_userparams(self):
        ffpar = {
            "RhoSq" : 1.207,
            "h_A1"  : 0.908,
            "R1"    : 1.401,
            "R2"    : 0.854,
            "R0"    : 1.15 ,
            "qiqj"  : "bc"   # Remember to override for b->u transitions!
        }
        return ffpar

    def isgur_wise(self, q2: float, ff: dict) -> dict:
        """Isgur-Wise-like relations to obtain tensor
        form factors from vector form factors from the equations of motion
        in the heavy quark limit

        Implements eq. 11 in https://arxiv.org/pdf/1503.05534 and assumes
        r(q2) = 1 for all three terms in the equation 

        Directly lifted from flavio https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/isgurwise.py

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at
        ff : dict
            dictionary with other FFs

        Returns
        -------
        dict
            FF dictionary with a tensor FFs
        """
        mB = self.get_param(f"m_{self.B}")
        mV = self.get_param(f"m_{self.V}")
        mb = running.get_mb_pole(self.pm)
        mq = running.get_mc_pole(self.pm) if self.get_userparam("qiqj") == "bc" else 0.
        # # power corrections
        # # NOTE: These are all zero
        # a_T1  = self.par[self._process + ' IW a_T1']
        # a_T2  = self.par[self._process + ' IW a_T2']
        # a_T23 = self.par[self._process + ' IW a_T23']
        a_T1, a_T2, a_T23 = 0.0, 0.0, 0.0
        # cf. eq. (11) of arXiv:1503.05534
        ff['T1'] = (mb + mq)/(mB + mV)*ff['V']  * ( 1 + a_T1 )
        ff['T2'] = (mb - mq)/(mB - mV)*ff['A1'] * ( 1 + a_T2 )
        if q2 == 0:
            ff['T23'] = (4*ff['A12']*(mb - mq)*(mB**2 + mV**2))/((mB - mV)**2*(mB + mV))
        else:
            ff['T23'] = ((mb - mq)* (
                    ((mB - mV)**2 - q2)* ((mB + mV)**2 - q2)*ff['A0']
                    + 8*mB*mV* (-mB**2 + mV**2)* ff['A12']
                    ))/ (4.*mB*(mV - mB)*mV*q2) * ( 1 + a_T23 )
        return ff

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Central value of $B\to V$ form factors in the lattice convention
        CLN parametrization. See eqs. (B4)-(B7) of arXiv:1203.2654.

        Directly lifted from flavio, see https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py

        For SM, this behaviour should be the same as HQET2 in EvtGen and CLN in HAMMER (though haven't verified throughly)
        
        For BSM, this follows flavio behaviour implementing eq. 11 in https://arxiv.org/pdf/1503.05534 and assuming
        r(q2) = 1 for all three terms in the equation 

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict[str, float]
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
        ff = self.isgur_wise(q2=q2, ff=ff)
        return ff