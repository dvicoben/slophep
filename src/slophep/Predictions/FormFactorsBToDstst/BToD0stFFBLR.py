import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD0st
from flavio.physics.bdecays.formfactors import hqet

class BLR_BToD0st(FormFactorBToD0st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD0st_BLR"
        self._ffpar = {
            "zt1"   : 0.7,
            "ztp"   : 0.2,
            "zeta1" : 0.6,
            "chi1"  : 0.0,
            "chi2"  : 0.0,
        }
        self._params = [
            "zt1", "ztp", "zeta1", "chi1", "chi2"
        ]
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD0starBLR.cc
        internalparams = {
            "ash" : 0.26/np.pi   ,
            "mb"  : 4.710        ,
            "mc"  : 4.710 - 3.400,
            "laB" : 0.4          ,
            "laS" : 0.76         ,
        }
        self._internalparams.update(internalparams)
    
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD0starBLR.cc
        """
        Mc = mC
        Mb = mB if mB else self.internalparams["Mb"]

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)
        mc = self.internalparams["mc"]
        mb = self.internalparams["mb"]
        zBC = mc/mb
        eB = 1./(2*mb)
        eC = 1./(2*mc)
        ash = self.internalparams["ash"]
        laB = self.internalparams["laB"]
        laS = self.internalparams["laS"]

        zt1   = self.ffpar["zt1"]
        ztp   = self.ffpar["ztp"]
        zeta1 = self.ffpar["zeta1"]
        chi1  = self.ffpar["chi1"]
        chi2  = self.ffpar["chi2"]

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        LOIWzeta = zt1 + (w-1)*zt1*ztp

        # QCD correction functions
        Cps = hqet.CP(w, zBC)
        # Cv1 = hqet.CV1(w, zBC)
        # Cv2 = hqet.CV2(w, zBC)
        # Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        # Ct2 = hqet.CT2(w, zBC)
        # Ct3 = hqet.CT3(w, zBC)
        
        gps = eC*(3*LambdaD12 - 2*(-1 + w*w)*zeta1 + (w-1)*(6*chi1 - 2*(1 + w)*chi2)) + (w-1)*(1 + ash*Cps) - (1 + w)*eB*Gb
        gp = -(eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1)) + ((w-1)*ash*(Ca2 + Ca3))/2. - eB*Gb
        gm = 1 + eC*(6*chi1 - 2*(1 + w)*chi2) + ash*(Ca1 + ((w-1)*(Ca2 - Ca3))/2.)
        gt = 1 + eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1 + 6*chi1 - 2*(1 + w)*chi2) + ash*Ct1 - eB*Gb

        ff = {
            "gP" : LOIWzeta*gps,
            "g+" : LOIWzeta*gp,
            "g-" : LOIWzeta*gm,
            "gT" : LOIWzeta*gt,
        }
        return ff

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])