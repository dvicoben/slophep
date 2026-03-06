import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1st
from flavio.physics.bdecays.formfactors import hqet

class BLR_BToD1st(FormFactorBToD1st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD1st_BLR"
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
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1starBLR.cc
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
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1starBLR.cc
        """
        Mc = mC
        Mb = mB if mB else self.internalparams["Mb"]

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)
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
        # Cs
        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        # Ca2 = hqet.CA2(w, zBC)
        # Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        # Gs  = 1 - eC*(LambdaD12/(1 + w) - 2*(w-1)*zeta1 + 2*chi1 - 2*(1 + w)*chi2) + as*Cs - eB*Gb;
        Gv1 = eC*(LambdaD12 - 2*(w-1)*chi1) + (w-1)*(1 + ash*Cv1) - (1 + w)*eB*Gb
        Gv2 = eC*(2*zeta1 - 2*chi2) - ash*Cv2
        Gv3 = -1 - eC*(LambdaD12/(1 + w) + 2*zeta1 - 2*chi1 + 2*chi2) - ash*(Cv1 + Cv3) + eB*Gb
        Ga  = 1 + eC*(LambdaD12/(1 + w) - 2*chi1) + ash*Ca1 - eB*Gb
        Gt1 = -1 + eC*(LambdaD12/(1 + w) + 2*chi1) - ash*(Ct1 + (w-1)*Ct2) + eB*Gb
        Gt2 = 1 + eC*(LambdaD12/(1 + w) - 2*chi1) + ash*(Ct1 - (w-1)*Ct3) + eB*Gb
        Gt3 = eC*(2*zeta1 + 2*chi2) - ash*Ct2

        ff = {
            "gS"  : 0.0         , #Gs, Need to fix for NP contributions 
            "gV1" : LOIWzeta*Gv1,
            "gV2" : LOIWzeta*Gv2,
            "gV3" : LOIWzeta*Gv3,
            "gA"  : LOIWzeta*Ga ,
            "gT1" : LOIWzeta*Gt1,
            "gT2" : LOIWzeta*Gt2,
            "gT3" : LOIWzeta*Gt3,
        }
        return ff

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])