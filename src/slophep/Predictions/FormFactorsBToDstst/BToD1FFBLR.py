import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1
from flavio.physics.bdecays.formfactors import hqet

class BLR_BToD1(FormFactorBToD1):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD1_BLR"
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
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
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
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1BLR.cc
        """
        Mc = mC
        Mb = mB if mB else self.internalparams["Mb"]

        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)
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

        LambdaD32 = -laB + laP*w
        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        LOIWtau = t1 + (w-1)*t1*tp

        # QCD correction functions
        Cv1 = hqet.CV1(w, zBC)
        Cv2 = hqet.CV2(w, zBC)
        Cv3 = hqet.CV3(w, zBC)
        Ca1 = hqet.CA1(w, zBC)
        # Ca2 = hqet.CA2(w, zBC)
        # Ca3 = hqet.CA3(w, zBC)
        Ct1 = hqet.CT1(w, zBC)
        Ct2 = hqet.CT2(w, zBC)
        Ct3 = hqet.CT3(w, zBC)

        # Fsc = (-(eC*(4*LambdaD32 + 2*(1 + w)*(6*eta1 + 2*(w-1)*eta2 - eta3) - 2*(w-1)*((1 + 2*w)*tau1 + tau2))) - 2*(1 + w)*(1 + ash*Cs) - 2*(w-1)*eB*Fb)/np.sqrt(6.)
        Fv1 = (-(eC*(4*(1 + w)*LambdaD32 - (-1 + w*w)*(2*eta1 + 3*eta3 + 3*tau1 - 3*tau2))) + (1 - w*w)*(1 + ash*Cv1) - (-1 + w*w)*eB*Fb)/np.sqrt(6.)
        Fv2 = (-3 - eC*(10*eta1 + 4*(w-1)*eta2 - 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) - ash*(3*Cv1 + 2*(1 + w)*Cv2) - 3*eB*Fb)/np.sqrt(6.)
        Fv3 = (-2 + w + eC*(4*LambdaD32 - 2*(6 + w)*eta1 - 4*(w-1)*eta2 - (-2 + 3*w)*eta3 + (2 + w)*tau1 + (2 + 3*w)*tau2) - ash*((2 - w)*Cv1 + 2*(1 + w)*Cv3) + (2 + w)*eB*Fb)/np.sqrt(6.)
        Fa = (-(eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*Ca1) - (w-1)*eB*Fb)/np.sqrt(6.)
        Ft1 = (-(eC*(4*LambdaD32 + (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (1 + w)*(1 + ash*(Ct1 + (w-1)*Ct2)) + (w-1)*eB*Fb)/np.sqrt(6.)
        Ft2 = (-(eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2))) + (-1 - w)*(1 + ash*(Ct1 - (w-1)*Ct3)) + (w-1)*eB*Fb)/np.sqrt(6.)
        Ft3 = (3 - eC*(-10*eta1 - 4*(w-1)*eta2 + 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) + ash*(3*Ct1 - (2 - w)*Ct2 + 3*Ct3) + 3*eB*Fb)/np.sqrt(6.)

        ff = {
            "fS"  : 0.0        , #Need to fix for NP contributions 
            "fV1" : LOIWtau*Fv1,
            "fV2" : LOIWtau*Fv2,
            "fV3" : LOIWtau*Fv3,
            "fA"  : LOIWtau*Fa ,
            "fT1" : LOIWtau*Ft1,
            "fT2" : LOIWtau*Ft2,
            "fT3" : LOIWtau*Ft3,
        }
        return ff

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])