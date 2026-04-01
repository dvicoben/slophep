import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD2st


class LLSW_BToD2st(FormFactorBToD2st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD2st_LLSW"
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
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD2starLLSW.cc
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

        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        # LOIWtau = t1 + tp*(w-1)  # Hammer
        LOIWtau = t1*(1. + tp*(w-1)) # basf2 EvtGenLLSWFF

        Ka1 = -1. - w - eC*(-((1 + w)*(2*eta1 - eta3)) + (w-1)*(tau1 - tau2)) - (w-1)*eB*Fb
        Ka2 = -2.*eC*(eta2 + tau1)
        Ka3 = 1. - eC*(2*eta1 - 2*eta2 - eta3 + tau1 + tau2) + eB*Fb
        Kv  = -1. - eC*(-2*eta1 + eta3 + tau1 - tau2) - eB*Fb

        ffs = {
            "kP"  : 0.0,
            "kA1" : LOIWtau*Ka1,
            "kA2" : LOIWtau*Ka2,
            "kA3" : LOIWtau*Ka3,
            "kV"  : LOIWtau*Kv ,
            "kT1" : 0.0,
            "kT2" : 0.0,
            "kT3" : 0.0
        }

        return ffs

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])