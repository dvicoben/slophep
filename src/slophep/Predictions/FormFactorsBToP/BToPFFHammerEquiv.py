from slophep.Predictions import FormFactorsBToP as FFBToP
import numpy as np

class BLPR_BToP_Hammer(FFBToP.BLPR_BToP):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale, *ffargs)
        self._name = "BToP_BLPR_Hammer"
        self._ffpar = {
            "RhoSq" : 1.24,
            "Chi21" : -0.06,
            "Chi2p" : 0.0,
            "Chi3p" : 0.05,
            "Eta1"  : 0.30,
            "Etap"  : -0.05,
            "dV20"  : 0.0
        }
        print(f"WARNING: Hammer uses a different basis for fT than the one assumed by flavio to compute observables. "
              +"Take care if using this for anything other than FF comparisons!")
    
    def get_ff(self, q2: float) -> dict:
        """FF in BLPR parameterisation from https://arxiv.org/pdf/1703.05330 as in HAMMER v1.2.1, 
        https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBLPR.cc?ref_type=tags

        Note that this returns fT=hT/sqrt(Mb*Mp) but this is not the basis flavio expects to compute observables.
        Take care if using for observables.

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict
            FF dictionary
        """
        ff = super().get_ff(q2)
        ff["fT"] = ff["fT"]/(self.internalparams["Mb"] + self.internalparams["Mc"])
        return ff


class BGL_BToP_Hammer(FFBToP.BGL_BToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)
        self._name = "BToP_BGL_Hammer"
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
        internalparams = {
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        }
        self._internalparams.update(internalparams)
    
    def get_ff(self, q2: float) -> dict:
        ff = super().get_ff(q2)
        ff["fT"] = ff["fT"]/(self.internalparams["Mb"] + self.internalparams["Mc"])
        return ff


class CLN_BToP_Hammer(FFBToP.CLN_BToP):
    def __init__(self, B: str, P: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, P, par, scale)
        self._name = "BToP_CLN_Hammer"
        # Parameters from B->D CLN Hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDCLN.cc
        self._ffpar = {
            "RhoSq" : 1.186,
            "G1"    : 1.082,
            "Delta" : 1.0
        }
        internalparams = {
            "a"      : 1.0,  # zero recoil expansion
            "ash"    : 0.26/np.pi,
            "mcOnmb" : 0.29,
            "lamBar" : 0.48
        }
        self._internalparams.update(internalparams)
    
    def get_ff(self, q2: float) -> dict:
        ff = super().get_ff(q2)
        ff["fT"] = ff["fT"]/(self.internalparams["Mb"] + self.internalparams["Mc"])
        return ff