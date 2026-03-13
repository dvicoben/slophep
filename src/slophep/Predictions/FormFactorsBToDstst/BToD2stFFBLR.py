import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD2st
from flavio.physics.bdecays.formfactors import hqet

class BLR_BToD2st(FormFactorBToD2st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD2st_BLR"
        self._ffpar = {
            "t1"   : 0.7 ,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.  ,
            "eta2" : 0.  ,
            "eta3" : 0.  ,
        }
        self._params = [
            "t1", "tp", "tau1", "tau2",
            "eta1", "eta2", "eta3"
        ]
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD2starBLR.cc
        internalparams = {
            "ash"  : 0.26/np.pi   ,
            "mb"   : 4.710        ,
            "mc"   : 4.710 - 3.400,
            "laB"  : 0.4          ,
            "laP"  : 0.8          ,
        }
        self._internalparams.update(internalparams)
    
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        """
        BLR form factors as in https://arxiv.org/pdf/1711.03110,
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD2starBLR.cc
        """
        Mc = mC
        Mb = mB if mB else self.internalparams["Mb"]

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = max((Mb2 + Mc2 - q2)/(2.*Mb*Mc), 1.0)
        zBC = self.internalparams["mc"]/self.internalparams["mb"]
        eB = 1./(2.*self.internalparams["mb"])
        eC = 1./(2.*self.internalparams["mc"])
        ash = self.internalparams["ash"]
        laB = self.internalparams["laB"]
        laP = self.internalparams["laP"]

        t1   = self.ffpar["t1"]
        tp   = self.ffpar["tp"]
        tau1 = self.ffpar["tau1"]
        tau2 = self.ffpar["tau2"]
        eta1 = self.ffpar["eta1"]
        eta2 = self.ffpar["eta2"]
        eta3 = self.ffpar["eta3"]

        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        LOIWtau = t1 + (w-1)*t1*tp

        # QCD correction functions
        # Cps = hqet.CP(w, zBC)
        Cv1 = hqet.CV1(w, zBC)
        # Cv2 = hqet.CV2(w, zBC)
        # Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        Ca2 = hqet.CA2(w, zBC)
        Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        # Kp = 1 + eC*(-2*eta1 - 2*(w-1)*eta2 + eta3 + (1 + 2*w)*tau1 + tau2) + as*Cps + eB*Fb
        Ka1 = -(eC*(-((1 + w)*(2*eta1 - eta3)) + (w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*Ca1) - (w-1)*eB*Fb
        Ka2 = -2*eC*(eta2 + tau1) + ash*Ca2
        Ka3 = 1 - eC*(2*eta1 - 2*eta2 - eta3 + tau1 + tau2) + ash*(Ca1 + Ca3) + eB*Fb
        Kv = -1 - eC*(-2*eta1 + eta3 + tau1 - tau2) - ash*Cv1 - eB*Fb
        Kt1 = 1 - eC*(2*eta1 - eta3) + ash*(Ct1 + ((w-1)*(Ct2 - Ct3))/2.)
        Kt2 = -(eC*(tau1 - tau2)) + ((1 + w)*ash*(Ct2 + Ct3))/2. + eB*Fb
        Kt3 = 2*eC*(-eta2 + tau1) - ash*Ct2

        ff = {
            "kP"  : 0.0        , #Need to fix for NP contributions 
            "kA1" : LOIWtau*Ka1,
            "kA2" : LOIWtau*Ka2,
            "kA3" : LOIWtau*Ka3,
            "kV"  : LOIWtau*Kv ,
            "kT1" : LOIWtau*Kt1,
            "kT2" : LOIWtau*Kt2,
            "kT3" : LOIWtau*Kt3,
        }
        return ff

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])