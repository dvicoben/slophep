from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToP import FormFactorBToP

class BGL_BToP(FormFactorBToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)

        self._name = "BToP_BGL"
        # Parameters from B->D BGL Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBGL.cc
        self._ffpar = {
            "f+_0" : 0.01565,
            "f+_1" : -0.0353,
            "f+_2" : -0.043,
            "f+_3" : 0.194,
            "f0_0" : 0.07932,
            "f0_1" : -0.214,
            "f0_2" : 0.17,
            "f0_3" : -0.958
        }
        self._params = [
            "f+_0", "f+_1", "f+_2", "f+_3",
            "f0_0", "f0_1", "f0_2", "f0_3"
        ]
        internalparams = {
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} Tensor FFs are 0, SM only parameterisation")

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr)
    
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

        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]
        Mb2 = Mb*Mb
        Mb3 = Mb2*Mb
        rC = Mc/Mb
        rC2 = rC*rC
        sqrC = np.sqrt(rC)

        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1)
        z = (sqrt(w+1) - sqrt(2))/(sqrt(w+1) + sqrt(2))
        zpow = np.array([z**ik for ik in range(int(self.internalparams["nmax"]))])

        ap = np.array([self.ffpar["f+_0"], self.ffpar["f+_1"], self.ffpar["f+_2"], self.ffpar["f+_3"]])
        a0 = np.array([self.ffpar["f0_0"], self.ffpar["f0_1"], self.ffpar["f0_2"], self.ffpar["f0_3"]])
        
        chiT = self.internalparams["ChiT"]
        chiL = self.internalparams["ChiL"]
        nc = self.internalparams["nc"]

        Pp = self.blaschke(self.internalparams["BcStatesp"], z, Mb, Mc)
        P0 = self.blaschke(self.internalparams["BcStates0"], z, Mb, Mc)

        phip = 32.*sqrt(nc/(6.*np.pi*chiT*Mb2))*rC2*pow(1+z,2)*pow(1-z,0.5)/pow((1+rC)*(1-z)+2*sqrC*(1+z), 5)
        phi0 = (1-rC2)*8.*sqrt(nc/(8.*np.pi*chiL))*rC*(1+z)*pow(1-z,1.5)/pow((1+rC)*(1-z)+2*sqrC*(1+z), 4)

        fp = np.dot(ap, zpow)/(Pp*phip)
        f0 = np.dot(a0, zpow)/(P0*phi0)

        ff = {
            "f+" : fp,
            "f0" : f0,
            "fT" : 0.0
        }
        return ff