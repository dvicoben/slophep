from slophep.Predictions import FormFactorsBToV as FFBToV
import numpy as np

class BGL_BToV_Hammer(FFBToV.BGL_BToV):
    def __init__(self, B: str, V: str, par: dict, scale: float, *ffargs):
        super().__init__(B, V, par, scale, *ffargs)
        self._name = "BToV_BGL_Hammer"
        self._ffpar = {
            "a0" : 0.00038,
            "a1" : 0.026905,
            "a2" : 0.,
            "b0" : 0.00055,
            "b1" : -0.0020370,
            "b2" : 0.,
            "c1" : -0.000433,
            "c2" : 0.005353,
            "d0" : 0.007,
            "d1" : -0.036
        }
        internalparams = {
            "nmax"       : 4,
            "Vcb"        : 41.5e-3,                       
            "chim"       : 3.068e-4, # GeV^-2
            "chip"       : 5.280e-4, # GeV^-2
            "chimL"      : 2.466e-3,
            "nc"         : 2.6,
            "etaEW"      : 1.0066,
            "BcStatesf"  : np.array([6.730, 6.736, 7.135, 7.142]),  # GeV
            "BcStatesg"  : np.array([6.337, 6.899, 7.012, 7.280]),  # GeV
            "BcStatesP1" : np.array([6.275, 6.842, 7.250])          # GeV
        }
        self._internalparams.update(internalparams)

        print("WARNING: Hammer's BToV BGL has a 1/(etaEW*Vcb) scaling of FFs, with etaEW=1.0066, Vcb=41.5e-3 and this results in different coefficients")
    
    def get_ff(self, q2: float) -> dict:
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags
        
        Note that there is an additional 1./(etaEW*Vcb) factor applied to FFs as in Hammer.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        etaEWVcb = self.internalparams["etaEW"]*self.internalparams["Vcb"]
        ff = super().get_ff(q2)
        return {k : ff[k]/etaEWVcb for k in ff}



class BLPR_BToV_Hammer(FFBToV.BLPR_BToV):
    def __init__(self, B: str, V: str, par: dict, scale: float, *ffargs):
        super().__init__(B, V, par, scale, *ffargs)
        self._name = "BToV_BLPR_Hammer"
        self._ffpar = {
            "RhoSq" : 1.24,
            "Chi21" : -0.06,
            "Chi2p" : 0.0,
            "Chi3p" : 0.05,
            "Eta1"  : 0.30,
            "Etap"  : -0.05,
            "dV20"  : 0.0
        }
        internalparams = {
            "ash"       : 0.26/np.pi,
            "la"        : 0.57115,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
            "rD"        : self.par['m_D0']/self.par['m_B0']
        }
        self._internalparams.update(internalparams)



class CLN_BToV_Hammer(FFBToV.CLN_BToV):
    def __init__(self, B: str, V: str, par: dict, scale: float, *ffargs):
        super().__init__(B, V, par, scale, *ffargs)
        self._name = "BToV_CLN_Hammer"
        self._ffpar = {
            "RhoSq" : 1.207,
            "h_A1"  : 0.908,
            "R1"    : 1.401,
            "R2"    : 0.854,
            "R0"    : 1.15
        }