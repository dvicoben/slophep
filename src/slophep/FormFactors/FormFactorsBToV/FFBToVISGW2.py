# Copyright (C) 2026  David Vico Benet

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# SLOP or SLOPHEP employs, translates and/or reimplements utilities from:
# - flavio (https://flav-io.github.io/), which is distributed under the MIT License, 
# and without any warranty, see <https://mit-license.org/>
# - Hammer (https://hammer.physics.lbl.gov/), which is distributed under version 3 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>
# - EOS (https://eoshep.org/), which is distributed under version 2 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>

import numpy as np
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)

class FFBToV_ISGW2(FormFactorBToV):
    """ISGW2 form-factors
    """
    _name = "FFBToV@ISGW2"
    
    def __init__(self, B: str, V: str):
        """ISGW2 form-factors
        """
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, V)

    def define_userparams(self):
        ffpar = {
            "msb" : 5.2,
            "msd" : 0.33,
            "bb2" : 0.431*0.431,
            "mbb" : 5.31,
            "nf"  : 4.0,
            "cf"  : 0.989,
            "msq" : 1.82,
            "bx2" : 0.38*0.38,
            "mbx" : 0.75*2.01+0.25*1.87,
            "nfp" : 3.0,
            "mqm" : 0.1
        }
        return ffpar

    def get_gammaji(self, z: float) -> float:
        value = 2+((2.0*z)/(1-z))*np.log(z)
        return -1.0*value

    def get_as(self, mq1: float, mq2: float) -> float:
        if mq2 <= 0.6:
            return 0.6
        lambdaSq = 0.04
        Nf = 4 if mq1 >= 1.85 else 3.0
        value = 12.0*np.pi/(33.0-2.0*Nf)
        value /= np.log(mq2*mq2/lambdaSq)
        return value

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """ISGW2 FFs

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict[str, float]
            FF dictionary
        """
        msb = self.get_userparam("msb")
        msd = self.get_userparam("msd")
        msq = self.get_userparam("msq")
        mbb = self.get_userparam("mbb")
        mbx = self.get_userparam("mbx")
        bb2 = self.get_userparam("bb2")
        bx2 = self.get_userparam("bx2")

        mtb=msb+msd
        mtx=msq+msd
        mup=1.0/(1.0/msq+1.0/msb)
        mum=1.0/(1.0/msq-1.0/msb)
        bbx2=0.5*(bb2+bx2)
        mb = self.get_param(f"m_{self.B}")
        mx = self.get_param(f"m_{self.V}")
        tm=(mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt=1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.get_userparam("mqm")
        nfp = self.get_userparam("nfp")

        r2 = 3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + (16.0/(mbb*mbx*(33.0-2.0*nfp)))*np.log(self.get_as(mqm,mqm) / self.get_as(msq,msq))
        ai = -1.0* ( 6.0/( 33.0 - 2.0*self.get_userparam("nf")))
        cji = np.power((self.get_as(msb,msb) / self.get_as(msq,msq)), ai)
        zji = msq / msb
        gammaji = self.get_gammaji(zji)
        chiji = -1.0 - ( gammaji / ( 1- zji ))
        
        betaji_g = (2.0/3.0)+gammaji
        betaji_f = (-2.0/3.0)+gammaji
        betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+(2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))
        betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-(2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+gammaji

        r_g = cji*(1+(betaji_g*self.get_as(msq, np.sqrt(mb*msq))/(np.pi)))
        r_f = cji*(1+(betaji_f*self.get_as(msq, np.sqrt(mb*msq))/(np.pi)))
        r_apmam = cji*(1+(betaji_apmam*self.get_as(msq, np.sqrt(mb*msq))/(np.pi)))
        
        f3  = np.sqrt(mtx/mtb)*np.power(np.sqrt(bx2*bb2)/bbx2, 1.5)/((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0))
        f3f = np.sqrt(mbx*mbb/(mtx*mtb))*f3
        f3g = np.sqrt(mtx*mtb/(mbx*mbb))*f3
        f3appam = np.sqrt(mtb*mtb*mtb*mbx/(mbb*mbb*mbb*mtx))*f3
        f3apmam = np.sqrt(mtx*mtb/(mbx*mbb))*f3

        cf = self.get_userparam("cf")
        f = cf*mtb*(1+wt+msd*(wt-1)/(2*mup))*f3f*r_f
        g = 0.5*(1/msq-msd*bb2/(2*mum*mtx*bbx2))*f3g*r_g

        appam = cji*(msd*bx2*(1-msd*bx2/(2*mtb*bbx2)) /
                     ((1+wt)*msq*msb*bbx2) - 
                     betaji_appam*self.get_as(msq, np.sqrt(msq*mb)) /
                     (mtb*np.pi))*f3appam
        apmam = -1.0*(mtb/msb-msd*bx2/(2*mup*bbx2)+wt*msd*mtb*bx2*
                      (1-msd*bx2/(2*mtb*bbx2))/((wt+1)*msq*msb*bbx2))*f3apmam*r_apmam/mtx
        ap = 0.5*(appam+apmam)
        am = 0.5*(appam-apmam)
        # MbcSqq = pow(mb + mx, 2.) - t
        # fs = -2 * mx * f / MbcSqq
        
        # We need to fiddle around with some conversions to go from hammer to flavio basis
        r  = mx/mb
        A1 = f/(mb + mx)
        pre = 1. / 2. / np.sqrt(mb * mx)
        A1_factor = pre * ((mb + mx)**2 - t) / (mb + mx)
        hA1 = A1/A1_factor
        hA2 = -(mb / np.sqrt(r))*(am + ap)
        hA3 = np.sqrt(r)*mb*(am - ap)
        A0 = pre * (((mb + mx)**2 - t) / (2 * mx) * hA1
                    - (mb**2 - mx**2 + t) / (2 * mb) * hA2
                    - (mb**2 - mx**2 - t) / (2 * mx) * hA3)
        A2 = pre * (mb + mx) * (hA3 + mx / mb * hA2)
        A12 = ((A1 * (mb + mx)**2 * (mb**2 - mx**2 - t)
                - A2 * (mb**4 + (mx**2 - t)**2
                - 2 * mb**2 * (mx**2 + t)))
                / (16. * mb * mx**2 * (mb + mx)))

        ff = {
            "V"   : (mb + mx) * g,
            "A0"  : A0,
            "A1"  : f/(mb + mx),
            "A2"  : A2,
            "A12" : A12,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "T3"  : 0.0,
            "T23" : 0.0
        }
        return ff