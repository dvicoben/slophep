from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToP import FormFactorBToP

class BCL_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)
        self._name = "BToP_BCL"
        # Parameters from B->Pi BCL Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoPiBCL.cc
        self._ffpar = {
            "f+_0" : 0.419 ,
            "f+_1" : -0.495,
            "f+_2" : -0.43 ,
            "f+_3" : 0.22  ,
            "f0_0" : 0.510 ,
            "f0_1" : -1.700,
            "f0_2" : 1.53  ,
            "f0_3" : 4.52  
        }
        self._params = [
            "f+_0", "f+_1", "f+_2", "f+_3",
            "f0_0", "f0_1", "f0_2", "f0_3"
        ]
        internalparams = {
            "m1m"    : 5.325,
            "q2cons" : False,  # impose f_+(q^2=0) = f_0(q^2=0)
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} Tensor FFs are 0, SM only parameterisation")


    def get_ff_mmeson(self, q2: float, mM: float, mB: float = None) -> dict[str, float]:
        # Follow calculation here: https://arxiv.org/pdf/1503.07839 - Eqn 5.3 and 5.4
        # and https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BCL/FFBtoPiBCL.cc
        Mb = self.internalparams["Mb"] if mB is None else mB
        Mu = mM

        m1m = self.internalparams["m1m"]
        P1m = 1-(q2/(m1m*m1m))

        # w = (Mb**2 + Mu**2 - q2) / (2. * Mb * Mu)
        tp = (Mb + Mu)*(Mb + Mu) if self.internalparams.get("tp") is None else self.internalparams.get("tp")
        tm = (Mb - Mu)*(Mb - Mu)
        t0 = tp*(1. - np.sqrt(1. - tm/tp))
        z = (np.sqrt(tp - q2) - np.sqrt(tp - t0))/(np.sqrt(tp - q2) + np.sqrt(tp - t0))
        z0 = (sqrt(tp) - sqrt(tp - t0))/(sqrt(tp) + sqrt(tp - t0))

        fp_vec = np.array([icoef for icoef in self.ffpar if "f+_" in icoef])
        fp_vec = np.array([self.ffpar[f"f+_{k}"] for k in range(len(fp_vec))])
        f0_vec = np.array([icoef for icoef in self.ffpar if "f0_" in icoef])
        f0_vec = np.array([self.ffpar[f"f0_{k}"] for k in range(len(f0_vec))])
        Nz = len(fp_vec)
        N0 = len(f0_vec)
        Nmax = max(Nz, N0)
        z_vec  = np.array([z**k for k in range(Nmax)])
        neg1_z = np.array([(-1.)**(n-Nz)*(n/Nz)*z_vec[Nz]] for n in range(Nz))

        f0 = np.sum(f0_vec*z_vec)
        fp = (1/P1m)*np.sum(
            fp_vec*(z_vec - neg1_z)
        )

        if self.internalparams.get("q2cons", False):
            z0_vec = np.array([z0**k for k in range(Nmax)])
            neg1_z0 = np.array([(-1.)**(n-Nz)*(n/Nz)*z0_vec[Nz]] for n in range(Nz))

            fpq2 = np.sum(fp_vec*(z_vec - neg1_z0))
            f0q2 = np.sum(f0_vec*z0_vec)
            f0 += (fpq2 - f0q2)/z0_vec[N0]*z_vec[N0]

        ff = {
            "f+" : fp ,
            "f0" : f0 ,
            "fT" : 0.0
        }
        return ff

    
    def get_ff(self, q2: float) -> dict:
        """Calculates BGL FFs.
        Implementation lifted from Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBGL.cc, 

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        return self.get_ff_mmeson(q2, self.internalparams["Mc"])