import numpy as np
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from slophep.FormFactors import hqet
from flavio.physics.bdecays.formfactors.b_v.cln import h_to_A

import logging
logger = logging.getLogger(__name__)

class FFBToV_BLPR(FormFactorBToV):
    _name = "FFBToV@BLPR"
    def define_userparams(self):
        ffpar = {
            "a"         : 1.509/np.sqrt(2),
            "rD"        : 1867./5280.,
            "RhoSq"     : 1.24      ,
            "Chi21"     : -0.06     ,
            "Chi2p"     : 0.0       ,
            "Chi3p"     : 0.05      ,
            "Eta1"      : 0.30      ,
            "Etap"      : -0.05     ,
            "dV20"      : 0.0       ,
            # internalparams:
            "mb"        : 4.710     ,
            "delta_mbc" : 3.4       ,
            "normscale" : 1.0       ,
            "ash"       : 0.26/np.pi,
            "ebReb"     : 0.861     ,
            "ecRec"     : 0.822     ,
            "mbarB"     : 5.313     ,
            "lam1"      : -0.3      , # GeV^2
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """FF in BLPR parameterisation from https://arxiv.org/pdf/1703.05330 as in HAMMER v1.4.1
        
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

        w = max(self.w(q2), 1)
        ash = self.get_userparam("ash")
        mb = self.get_userparam("mb")
        la = self.get_userparam("mbarB") - mb + self.get_userparam("lam1")/(2*mb)
        eb = la/(2*mb)
        mc = mb - self.get_userparam("delta_mbc")
        ec = la/(2*mc)
        ebReb = self.get_userparam("ebReb")
        ecRec = self.get_userparam("ecRec")
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)
        zBC   = mc/mb

        RhoSq = self.get_userparam("RhoSq")
        chi21 = self.get_userparam("Chi21")
        chi2p = self.get_userparam("Chi2p")
        chi3p = self.get_userparam("Chi3p")
        eta1  = self.get_userparam("Eta1" )
        etap  = self.get_userparam("Etap" )

        L1 = -4.0*(w-1.0)*(chi21 + (w-1.0)*chi2p)+12.0*chi3p*(w-1.0)
        L2 = -4.0*chi3p*(w-1.0)
        L3 = 4.0*(chi21 + (w-1.0)*chi2p)
        L4 = 2.*(eta1 + etap*(w-1.))-1.
        # L4 = (2.0*(eta1 + etap*(w-1.0))-1.0*ebReb)
        L5 = -1.0
        # L5 = -ecRec
        L6 = -2.*(1+(eta1 + etap*(w-1.)))/(w+1.)
        # L6 = -2.0*(ecRec + (eta1 + etap*(w-1.0)))/(w+1.0)

        # QCD correction functions
        Cv1 = hqet.CV1(w, zBC)
        # CV2 = hqet.CV2(w, zBC)
        # CV3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        # From Sec. III-A in https://arxiv.org/pdf/1703.05330 - Variables for leading IW function (derived from G(1))
        # rD = self.get_param('m_D0')/self.get_param('m_B0')
        # a = np.sqrt((1+rD)/(2*np.sqrt(rD)))
        # # Values differ slightly in hammer source, where a=1.509 / sqrt2 and rD=1867./5280.:
        rD = self.get_userparam("rD")
        a = self.get_userparam("a")
        V21 = 57.0
        V20 = 7.5 + self.get_userparam("dV20")
        zCon = (np.sqrt(w+1.0) - np.sqrt(2.0)*a)/(np.sqrt(w+1.0) + np.sqrt(2.0)*a)

        # Derivatives (approx) at w_0 = 2a^2-1
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

        # Hps = 1.+ash*Cps+ec*(L2+L3*(w-1)+L5-L6*(w+1))+eb*(L1-L4)-(corrc + corrb)
        Hv = 1.+ash*Cv1+ec*(L2-L5)+eb*(L1-L4) -(corrc + corrb)
        Ha1 = 1.+ash*Ca1+ec*(L2-L5*(w-1)/(w+1))+eb*(L1-L4*(w-1)/(w+1)) -(corrc + corrb)*(w-1)/(w+1)
        Ha2 = ash*Ca2+ec*(L3+L6) +2.*corrc/(w+1.)
        Ha3 = 1.+ash*(Ca1+Ca3)+ec*(L2-L3-L5+L6)+eb*(L1-L4) -(corrc * (w-1.)/(w+1.)+corrb)
        Ht1 = 1.+ash*(Ct1+0.5*(w-1)*(Ct2-Ct3))+ec*(L2)+eb*(L1)
        Ht2 = 0.5*(w+1)*ash*(Ct2+Ct3)+ec*(L5)-eb*(L4) +corrc-corrb
        Ht3 = ash*(Ct2)+ec*(L6-L3) +2.*corrc/(w+1.)
        
        normscale = self.get_userparam("normscale")
        h = {
            "V"  : normscale*xi*Hv,
            "A1" : normscale*xi*Ha1,
            "A2" : normscale*xi*Ha2,
            "A3" : normscale*xi*Ha3,
            "T1" : normscale*xi*Ht1,
            "T2" : normscale*xi*Ht2,
            "T3" : normscale*xi*Ht3
        }
        # NOTE: this performs the translation https://arxiv.org/pdf/1309.0301 eqns. 38-39,
        # should be analgous to eqns B7-B13 in https://arxiv.org/pdf/1908.09398
        return h_to_A(mB, mV, h, q2)
        # return ff
    
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
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)
        zBC   = mc/mb

        RhoSq = self.get_userparam("RhoSq")
        chi21 = self.get_userparam("Chi21")
        chi2p = self.get_userparam("Chi2p")
        chi3p = self.get_userparam("Chi3p")
        eta1  = self.get_userparam("Eta1" )
        etap  = self.get_userparam("Etap" )

        L1 = -4.0*(w-1.0)*(chi21 + (w-1.0)*chi2p)+12.0*chi3p*(w-1.0)
        L2 = -4.0*chi3p*(w-1.0)
        L3 = 4.0*(chi21 + (w-1.0)*chi2p)
        L4 = 2.*(eta1 + etap*(w-1.))-1.
        # L4 = (2.0*(eta1 + etap*(w-1.0))-1.0*ebReb)
        L5 = -1.0
        # L5 = -ecRec
        L6 = -2.*(1+(eta1 + etap*(w-1.)))/(w+1.)
        # L6 = -2.0*(ecRec + (eta1 + etap*(w-1.0)))/(w+1.0)

        # QCD correction functions
        Cv1 = hqet.CV1(w, zBC)
        # CV2 = hqet.CV2(w, zBC)
        # CV3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)
        
        # Hps = 1.+ash*Cps+ec*(L2+L3*(w-1)+L5-L6*(w+1))+eb*(L1-L4)-(corrc + corrb)
        Hv = 1.+ash*Cv1+ec*(L2-L5)+eb*(L1-L4) -(corrc + corrb)
        Ha1 = 1.+ash*Ca1+ec*(L2-L5*(w-1)/(w+1))+eb*(L1-L4*(w-1)/(w+1)) -(corrc + corrb)*(w-1)/(w+1)
        Ha2 = ash*Ca2+ec*(L3+L6) +2.*corrc/(w+1.)
        Ha3 = 1.+ash*(Ca1+Ca3)+ec*(L2-L3-L5+L6)+eb*(L1-L4) -(corrc * (w-1.)/(w+1.)+corrb)
        Ht1 = 1.+ash*(Ct1+0.5*(w-1)*(Ct2-Ct3))+ec*(L2)+eb*(L1)
        Ht2 = 0.5*(w+1)*ash*(Ct2+Ct3)+ec*(L5)-eb*(L4) +corrc-corrb
        Ht3 = ash*(Ct2)+ec*(L6-L3) +2.*corrc/(w+1.)
        
        hhat = {
            "V"  : Hv,
            "A1" : Ha1,
            "A2" : Ha2,
            "A3" : Ha3,
            "T1" : Ht1,
            "T2" : Ht2,
            "T3" : Ht3
        }
        return hhat