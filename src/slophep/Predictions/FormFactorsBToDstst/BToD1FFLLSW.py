import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1


class LLSW_BToD1(FormFactorBToD1):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD1_LLSW"
        self._ffpar = {
            "t1"   : 0.71,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.0 ,
            "eta2" : 0.0 ,
            "eta3" : 0.0 ,
        }
        self._params = []
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1LLSW.cc
        internalparams = {
            "mb"  : 4.2,
            "mc"  : 1.4,
            "laB" : 0.4,
            "laP" : 0.8 
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} SM only, tensor FFs are zero.")

    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = self.internalparams["Mb"] if mB is None else mB
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)

        eB   = 0.5 / self.internalparams["mb"]
        eC   = 0.5 / self.internalparams["mc"]
        laB  = self.internalparams["laB"]
        laP  = self.internalparams["laP"]
        t1   = self.ffpar["t1"]
        tp   = self.ffpar["tp"]
        tau1 = self.ffpar["tau1"]
        tau2 = self.ffpar["tau2"]
        eta1 = self.ffpar["eta1"]
        eta2 = self.ffpar["eta2"]
        eta3 = self.ffpar["eta3"]

        LambdaD32 = -laB + laP*w
        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        # LOIWtau = t1 + tp*(w-1)  # Hammer
        LOIWtau = t1*(1. + tp*(w-1)) # basf2 EvtGenLLSWFF

        Fv1 = (1 - w*w - eC*(4*(1 + w)*LambdaD32 - (-1 + w*w)*(2*eta1 + 3*eta3 + 3*tau1 - 3*tau2)) - (-1 + w*w)*eB*Fb)/np.sqrt(6.)
        Fv2 = (-3 - eC*(10*eta1 + 4*(w-1)*eta2 - 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) - 3*eB*Fb)/np.sqrt(6.)
        Fv3 = (-2 + w + eC*(4*LambdaD32 - 2*(6 + w)*eta1 - 4*(w-1)*eta2 - (-2 + 3*w)*eta3 + (2 + w)*tau1 + (2 + 3*w)*tau2) + (2 + w)*eB*Fb)/np.sqrt(6.)
        Fa = (-1 - w - eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2)) - (w-1)*eB*Fb)/np.sqrt(6.)
        
        ffs = {
            "fS"  : 0.0,
            "fV1" : LOIWtau*Fv1,
            "fV2" : LOIWtau*Fv2,
            "fV3" : LOIWtau*Fv3,
            "fA"  : LOIWtau*Fa ,
            "fT1" : 0.0,
            "fT2" : 0.0,
            "fT3" : 0.0
        }

        return ffs

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])