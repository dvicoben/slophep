from math import sqrt
import numpy as np

from flavio.physics.bdecays.formfactors import hqet
from flavio.physics.bdecays.formfactors.b_p.cln import h_to_f
from slophep.Predictions.FormFactorsBToP import FormFactorBToP

class BLPR_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)

        self._name = "BToP_BLPR"
        self._ffpar = {
            "RhoSq"     : 1.24,
            "Chi21"     : -0.06,
            "Chi2p"     : 0.0,
            "Chi3p"     : 0.05,
            "Eta1"      : 0.30,
            "Etap"      : -0.05,
            "dV20"      : 0.0,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "normscale" : 1.0
        }
        self._params = [
            "RhoSq", "Chi21", "Chi2p", "Chi3p", "Eta1", "Etap", "dV20",
            "mb", "delta_mbc", "normscale"
        ]
        internalparams = {
            "ash"       : 0.26/np.pi,
            "mbarB"     : 5.313,
            "lam1"      : -0.3, # GeV^2
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
            "rD"        : self.par['m_D0']/self.par['m_B0']
        }
        self._internalparams.update(internalparams)

    
    def get_ff(self, q2: float) -> dict:
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
        mB = self.internalparams["Mb"]
        mP = self.internalparams["Mc"]

        w = max(self.w(q2), 1)
        ash = self.internalparams["ash"]
        mb = self.ffpar["mb"]
        la = self.internalparams["mbarB"] - mb + self.internalparams["lam1"]/(2*mb)
        eb = la/(2*mb)
        mc = mb - self.ffpar["delta_mbc"]
        ec = la/(2*mc)
        ebReb = self.internalparams["ebReb"]
        ecRec = self.internalparams["ecRec"]
        zBC = mc/mb
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)

        RhoSq = self.ffpar["RhoSq"]
        chi21 = self.ffpar["Chi21"]
        chi2p = self.ffpar["Chi2p"]
        chi3p = self.ffpar["Chi3p"]
        eta1 = self.ffpar["Eta1"]
        etap = self.ffpar["Etap"]
        
        L1 = -4.*(w-1)*(chi21 + (w-1.)*chi2p)+12.*chi3p*(w-1.)
        L4 = 2.*(eta1 + etap*(w-1.))-1.
        # L4b = (2.*(eta1 + etap*(w-1.))-1.*ebReb) # HAMMERv1.2.1
        # L4c = (2.*(eta1 + etap*(w-1.))-1.*ecRec) # HAMMERv1.2.1

        rD = self.internalparams["rD"]
        a = sqrt((1+rD)/(2*sqrt(rD)))
        V21 = 57.0
        V20 = 7.5 + self.ffpar["dV20"]
        zCon = (np.sqrt(w+1.0) - sqrt(2.0)*a)/(np.sqrt(w+1.0) + sqrt(2.0)*a)

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

        normscale = self.ffpar["normscale"]
        ff = {
            "h+" : normscale*xi*Hp,
            "h-" : normscale*xi*Hm,
            "hT" : normscale*xi*Ht
        }
        return h_to_f(mB, mP, ff, q2)
    
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
        zBC = mc/mb
        corrb = eb*(1.-ebReb)
        corrc = ec*(1.-ecRec)

        chi21 = self.ffpar["Chi21"]
        chi2p = self.ffpar["Chi2p"]
        chi3p = self.ffpar["Chi3p"]
        eta1 = self.ffpar["Eta1"]
        etap = self.ffpar["Etap"]
        
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