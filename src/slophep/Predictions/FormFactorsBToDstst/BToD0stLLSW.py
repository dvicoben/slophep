import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD0st


class LLSW_BToD0st(FormFactorBToD0st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD0st_LLSW"
        self._ffpar = {
            "zt1"  : 0.68 ,
            "ztp"  : -0.2 ,
            "zeta1": 0.3  ,
        }
        self._params = []
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD1starLLSW.cc
        internalparams = {
            "mb"   : 4.2  ,
            "mc"   : 1.4  ,
            "laB"  : 0.4  ,
            "laS"  : 0.76 ,
            "chi1" : 0.   ,
            "chi2" : 0.   ,
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} SM only, tensor FFs are zero.")

    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = self.internalparams["Mb"] if mB is None else mB
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)
        eB = 0.5 / self.internalparams["mb"]
        eC = 0.5 / self.internalparams["mc"]
        zt1 = self.ffpar["zt1"]
        ztp = self.ffpar["ztp"]
        zeta1 = self.ffpar["zeta1"]
        chi1 = self.internalparams["chi1"]
        chi2 = self.internalparams["chi2"]
        laB = self.internalparams["laB"]
        laS = self.internalparams["laS"]

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        # The Belle II EvtGen has callculation of LO IW that differs from the EvtGen code
        # in https://github.com/belle2/basf2/blob/main/generators/evtgen/models/src/EvtLLSWFF.cc
        # LOIWzeta = zt1 + (w-1)*ztp
        LOIWzeta = zt1*(1. + (w-1)*ztp)

        gp = -(eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1)) - eB*Gb
        gm = 1 + eC*(6*chi1 - 2*(1 + w)*chi2)


        ffs = {
            "gP" : 0.0,
            "g+" : LOIWzeta*gp,
            "g-" : LOIWzeta*gm,
            "gT" : 0.0
        }

        return ffs

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])