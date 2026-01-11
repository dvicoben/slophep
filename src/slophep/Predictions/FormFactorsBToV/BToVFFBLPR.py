from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToV import FormFactorBToV
from flavio.physics.bdecays.formfactors import hqet
from flavio.physics.bdecays.formfactors.b_v.cln import h_to_A


class BLPR_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_BLPR"
        self._ffpar = {
            "RhoSq" : 1.24,
            "Chi21" : -0.06,
            "Chi2p" : 0.0,
            "Chi3p" : 0.05,
            "Eta1"  : 0.30,
            "Etap"  : -0.05,
            "dV20"  : 0.0
        }
        self._params = ["RhoSq", "Chi21", "Chi2p", "Chi3p", "Eta1", "Etap", "dV20"]
        internalparams = {
            "ash"       : 0.26/np.pi,
            "la"        : 0.57115,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
            "rD"        : self.par['m_D0']/self.par['m_B0']
        }
        self._internalparams.update(internalparams)

    def get_ff(self, q2: float) -> dict:
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
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]

        w = max(self.w(q2), 1)
        ash = self.internalparams["ash"]
        la = self.internalparams["la"]
        mb = self.internalparams["mb"]
        eb = la/(2*mb)
        mc = mb - self.internalparams["delta_mbc"]
        ec = la/(2*mc)
        ebReb = self.internalparams["ebReb"]
        ecRec = self.internalparams["ecRec"]
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)
        zBC = mc/mb

        # z = ffcommon.z(mB, mV, q2, t0='tm')
        RhoSq = self.ffpar["RhoSq"]
        chi21 = self.ffpar["Chi21"]
        chi2p = self.ffpar["Chi2p"]
        chi3p = self.ffpar["Chi3p"]
        eta1 = self.ffpar["Eta1"]
        etap = self.ffpar["Etap"]

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
        rD = self.internalparams["rD"]
        a = sqrt((1+rD)/(2*sqrt(rD)))
        V21 = 57.0
        V20 = 7.5 + self.ffpar["dV20"]
        zCon = (np.sqrt(w+1.0) - sqrt(2.0)*a)/(np.sqrt(w+1.0) + sqrt(2.0)*a)

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

        h = {
            "V" : xi*Hv,
            "A1" : xi*Ha1,
            "A2" : xi*Ha2,
            "A3" : xi*Ha3,
            "T1" : xi*Ht1,
            "T2" : xi*Ht2,
            "T3" : xi*Ht3
        }
        # NOTE: this performs the translation https://arxiv.org/pdf/1309.0301 eqns. 38-39,
        # should be analgous to eqns B7-B13 in https://arxiv.org/pdf/1908.09398
        return h_to_A(mB, mV, h, q2)

    def get_hhat(self, q2: float) -> dict:
        w = max(self.w(q2), 1)
        ash = self.internalparams["ash"]
        la = self.internalparams["la"]
        mb = self.internalparams["mb"]
        eb = la/(2*mb)
        mc = mb - self.internalparams["delta_mbc"]
        ec = la/(2*mc)
        ebReb = self.internalparams["ebReb"]
        ecRec = self.internalparams["ecRec"]
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)
        zBC = mc/mb

        chi21 = self.ffpar["Chi21"]
        chi2p = self.ffpar["Chi2p"]
        chi3p = self.ffpar["Chi3p"]
        eta1 = self.ffpar["Eta1"]
        etap = self.ffpar["Etap"]

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