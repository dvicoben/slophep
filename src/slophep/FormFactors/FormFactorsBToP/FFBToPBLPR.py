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
from slophep.FormFactors.FormFactorsBToP.FFBToPBase import FormFactorBToP
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from slophep.FormFactors import hqet
from flavio.physics.bdecays.formfactors.b_p.cln import h_to_f

import logging
logger = logging.getLogger(__name__)

class FFBToP_BLPR(FormFactorBToP):
    _name = "FFBToP@BLPR"
    def define_userparams(self):
        ffpar = {
            "a"         : 1.509/np.sqrt(2),
            "rD"        : 1867./5280.,
            "RhoSq"     : 1.24,
            "Chi21"     : -0.06,
            "Chi2p"     : 0.0,
            "Chi3p"     : 0.05,
            "Eta1"      : 0.30,
            "Etap"      : -0.05,
            "dV20"      : 0.0,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "normscale" : 1.0,
            "ash"       : 0.26/np.pi,
            "mbarB"     : 5.313,
            "lam1"      : -0.3 , # GeV^2
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """FF in BLPR parameterisation from https://arxiv.org/pdf/1703.05330 as in HAMMER v1.2.1, 
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBLPR.cc?ref_type=tags

        Note that HAMMER's basis differs from the one used here for fT, where in hammer
        fT=hT/sqrt(Mb*Mp) while here we use Appendix B in https://arxiv.org/abs/1908.09398, which differs
        by a factor of (Mb+Mp).

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
        mP = self.get_param(f"m_{self.P}")

        w = max(self.w(q2), 1)
        ash = self.get_userparam("ash")
        mb = self.get_userparam("mb")
        la = self.get_userparam("mbarB") - mb + self.get_userparam("lam1")/(2*mb)
        eb = la/(2*mb)
        mc = mb - self.get_userparam("delta_mbc")
        ec = la/(2*mc)
        ebReb = self.get_userparam("ebReb")
        ecRec = self.get_userparam("ecRec")
        zBC = mc/mb
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)

        RhoSq = self.get_userparam("RhoSq")
        chi21 = self.get_userparam("Chi21")
        chi2p = self.get_userparam("Chi2p")
        chi3p = self.get_userparam("Chi3p")
        eta1  = self.get_userparam("Eta1")
        etap  = self.get_userparam("Etap")
        
        L1 = -4.*(w-1)*(chi21 + (w-1.)*chi2p)+12.*chi3p*(w-1.)
        L4 = 2.*(eta1 + etap*(w-1.))-1.
        # L4b = (2.*(eta1 + etap*(w-1.))-1.*ebReb) # HAMMERv1.2.1
        # L4c = (2.*(eta1 + etap*(w-1.))-1.*ecRec) # HAMMERv1.2.1

        # # From Sec. III-A in https://arxiv.org/pdf/1703.05330 - Variables for leading IW function (derived from G(1))
        # rD = self.get_param("m_D0")/self.get_param("m_B0")
        # a = np.sqrt((1+rD)/(2*np.sqrt(rD)))
        # # Values differ slightly in hammer source, where a=1.509 / sqrt2 and rD=1867./5280.:
        rD = self.get_userparam("rD")
        a = self.get_userparam("a")
        V21 = 57.0
        V20 = 7.5 + self.get_userparam("dV20")
        zCon = (np.sqrt(w+1.0) - np.sqrt(2.0)*a)/(np.sqrt(w+1.0) + np.sqrt(2.0)*a)

        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        w0 = 2.0*a*a - 1.0
        Cv1z = hqet.CV1(w0, zBC)
        Cv2z = hqet.CV2(w0, zBC)
        Cv3z = hqet.CV3(w0, zBC)
        Cv1zp = (hqet.CV1(w0 + 1e-5, zBC) - Cv1z)/1e-5
        Cv2zp = (hqet.CV2(w0 + 1e-5, zBC) - Cv2z)/1e-5
        Cv3zp = (hqet.CV3(w0 + 1e-5, zBC) - Cv3z)/1e-5
        Cv1zpp = (Cv1zp - (Cv1z - hqet.CV1(w0 - 1e-5, zBC))/1e-5)/1e-5
        Cv2zpp = (Cv2zp - (Cv2z - hqet.CV2(w0 - 1e-5, zBC))/1e-5)/1e-5
        Cv3zpp = (Cv3zp - (Cv3z - hqet.CV3(w0 - 1e-5, zBC))/1e-5)/1e-5

        a2 = a*a
        a4 = a2*a2
        a6 = a4*a2
        Xi = 64*a4*RhoSq - 16*a2 - V21
        xiIW = 1 - 8*a2*RhoSq*zCon + zCon*zCon*(
                V21*RhoSq - V20 + (eb - ec)*(2*Xi*etap * (1-rD)/(1+rD))
                + (eb + ec)*(Xi* (12*chi3p - 4*chi21) - 16*((a2-1)*Xi - 16* a4)*chi2p)
                + ash*(Xi*(Cv1zp +(Cv3z + rD*Cv2z)/(1 + rD)) + 2*a2*(Xi - 32*a2)*(Cv3zp + rD*Cv2zp)/(1+rD) - 64*a6*(Cv3zpp + rD*Cv2zpp)/(1+rD) -32*a4*Cv1zpp ))
        xiIW /= 1 - 8*a2*RhoSq*(1-a)/(1+a) + pow((1-a)/(1+a),2.)*(
                V21*RhoSq - V20 + (eb - ec)*(2*Xi*etap * (1-rD)/(1+rD))
                + (eb + ec)*(Xi* (12*chi3p - 4*chi21) - 16*((a2-1)*Xi - 16* a4)*chi2p)
                + ash*(Xi*(Cv1zp +(Cv3z + rD*Cv2z)/(1 + rD)) + 2*a2*(Xi - 32*a2)*(Cv3zp + rD*Cv2zp)/(1+rD) - 64*a6*(Cv3zpp + rD*Cv2zpp)/(1+rD) -32*a4*Cv1zpp ))

        chi1=0
        xi =xiIW + 2.*(eb+ec)*chi1
        Hp = 1.+ash*(Cv1+0.5*(w+1)*(Cv2+Cv3))+(ec+eb)*L1
        # Hm = ash*0.5*(w+1)*(Cv2-Cv3)+(ec*L4c-eb*L4b)               # HAMMERv1.2.1
        Hm = ash*0.5*(w+1)*(Cv2-Cv3)+(ec-eb)*L4 + (corrc-corrb)
        # Ht = 1.+ash*(Ct1-Ct2+Ct3)+(ec+eb)*L1 - (ec*L4c+eb*L4b)     # HAMMERv1.2.1
        Ht = 1.+ash*(Ct1-Ct2+Ct3)+(ec+eb)*(L1 - L4) - (corrc+corrb)

        normscale = self.get_userparam("normscale")
        ff = {
            "h+" : normscale*xi*Hp,
            "h-" : normscale*xi*Hm,
            "hT" : normscale*xi*Ht
        }
        return h_to_f(mB, mP, ff, q2)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_hhat(self, q2: float) -> dict[str, float]:
        w = max(self.w(q2), 1)
        ash = self.get_userparam("ash")
        mb = self.get_userparam("mb")
        la = self.get_userparam("mbarB") - mb + self.get_userparam("lam1")/(2*mb)
        eb = la/(2*mb)
        mc = mb - self.get_userparam("delta_mbc")
        ec = la/(2*mc)
        ebReb = self.get_userparam("ebReb")
        ecRec = self.get_userparam("ecRec")
        zBC = mc/mb
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)

        RhoSq = self.get_userparam("RhoSq")
        chi21 = self.get_userparam("Chi21")
        chi2p = self.get_userparam("Chi2p")
        chi3p = self.get_userparam("Chi3p")
        eta1  = self.get_userparam("Eta1")
        etap  = self.get_userparam("Etap")
        
        L1 = -4.*(w-1)*(chi21 + (w-1.)*chi2p)+12.*chi3p*(w-1.)
        L4 = 2.*(eta1 + etap*(w-1.))-1.

        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        Hp = 1.+ash*(Cv1+0.5*(w+1)*(Cv2+Cv3))+(ec+eb)*L1
        Hm = ash*0.5*(w+1)*(Cv2-Cv3)+(ec-eb)*L4 + (corrc-corrb)
        Ht = 1.+ash*(Ct1-Ct2+Ct3)+(ec+eb)*(L1 - L4) - (corrc+corrb)

        h = {
            "h+" : Hp,
            "h-" : Hm,
            "hT" : Ht
        }
        return h