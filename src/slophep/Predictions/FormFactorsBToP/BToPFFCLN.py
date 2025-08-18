from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToP import FormFactorBToP

class CLN_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)

        self._name = "BToP_CLN"
        # Parameters from B->D CLN Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDCLN.cc
        self._ffpar = {
            "RhoSq" : 0.01565,
            "G1"    : -0.0353,
            "Delta" : 1.0
        }
        self._params = [
            "RhoSq", "G1"
        ]
        internalparams = {
            "a"      : 1.0,  # zero recoil expansion
            "ash"    : 0.26/np.pi,
            "mcOnmb" : 0.29,
            "lamBar" : 0.48
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} Tensor FFs are 0, SM only parameterisation")

    def get_ff(self, q2: float) -> dict:
        """Calculates CLN FFs.
        Implementation lifted from Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDCLN.cc, 

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """

        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]
        rC = Mc/Mb
        sqrt2 = sqrt(2)

        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1)
        a = self.internalparams["a"]
        zCon = (sqrt(w+1) - sqrt2*a)/(sqrt(w+1) + sqrt2*a)
        
        RhoSq = self.ffpar["RhoSq"]
        G1par = self.ffpar["G1"]
        Delta = self.ffpar["Delta"]

        V1wOnV1 = 1. - 8.* RhoSq * zCon + (51. * RhoSq - 10.) * zCon * zCon - (252. * RhoSq - 84.) * zCon * zCon *zCon
        fp=((1. + rC)/(2.*sqrt(rC))*V1wOnV1)*G1par
        S1wOnV1w = 1 + Delta*(-0.019 + 0.041*(w - 1.) - 0.015*pow(w - 1.,2.))
        f0 = (sqrt(rC)/(1. + rC)*(1+w)*S1wOnV1w*V1wOnV1)*G1par

        ff = {
            "f+" : fp,
            "f0" : f0,
            "fT" : 0.0
        }
        return ff