import numpy as np
from typing import Any
from slophep.FormFactors.FormFactorsBToV.FFBToVBase import FormFactorBToV
from slophep.Tools.SamplingTools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


def outer_function_vec(z: float, r: float, a: int, b: int, c: int) -> float:
    """missing a prefactor that is dependent on the exact FF"""
    phi = np.power(r*(1+z), a) * np.power(1-z, 0.5*b) / np.power(
        (1+r)*(1-z) + 2*np.sqrt(r)*(1+z), c
    )
    return phi


def outer_function_ten(q2: float, tp: float, tm: float, t0: float, a: int, b: int, c: int) -> float:
    """missing a prefactor that is dependent on the exact FF"""
    phi = (
        np.power((tp - q2)/(tp - t0), 0.25)
        *(np.sqrt(tp - q2) + np.sqrt(tp-t0))
        *np.power(tp - q2, 0.25*a)
        *np.power(np.sqrt(tp - q2) + np.sqrt(tp-tm), 0.5*b)
        *np.power(np.sqrt(tp - q2) + np.sqrt(tp), -(c+3))
    )
    return phi


class FFBToV_BGLGeneric(FormFactorBToV):
    _name = "FFBToV@BGLGen"
    
    def __init__(self, B: str, V: str,
                 N_g   : int = 3, 
                 N_f   : int = 3, 
                 N_F1  : int = 3, 
                 N_F2  : int = 3,
                 N_T1  : int = 0, 
                 N_T2  : int = 0, 
                 N_T23 : int = 0):
        """Generalised B->V BGL form factor computation

        Parameters
        ----------
        B : str
            B meson
        V : str
            Vector meson
        N_g : int, optional
            Number of coefficients in g, by default 3
        N_f : int, optional
            Number of coefficients in f, by default 3
        N_F1 : int, optional
            Number of coefficients in F1, by default 3
        N_F2 : int, optional
            Number of coefficients in F2, by default 2
        N_T1 : int, optional
            Number of coefficients in T1, by default 0
        N_T2 : int, optional
            Number of coefficients in T2, by default 0
        N_T23 : int, optional
            Number of coefficients in T23, by default 0
        """
        self._n = {
            "g"    : N_g  ,
            "f"    : N_f  ,
            "F1"   : N_F1 ,
            "F2"   : N_F2 ,
            "T1"   : N_T1 ,
            "T2"   : N_T2 ,
            "T23"  : N_T23,
            "nmax" : max(N_g, N_f, N_F1, N_F2, N_T1, N_T2, N_T23)
        }
        super().__init__(B, V)

    @property
    def n(self) -> dict[str, int]:
        return self._n

    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "chig"             : 5.131e-4 ,
            "chif"             : 3.894e-4 ,
            "chiF1"            : 3.894e-4 ,
            "chiF2"            : 1.9421e-2,
            "chiT1"            : 4.98e-4  ,
            "chiT2"            : 2.77e-4  ,
            "chiT23"           : 2.77e-4  ,
            "nc"               : 2.6      ,
            # "etaEW"            : 1.0066   ,
            "BcStates1p"       : np.array([6.739, 6.750, 7.145, 7.150]), #A1
            "BcStates1m"       : np.array([6.329, 6.910, 7.020]),        #V1
            "BcStates0m"       : np.array([6.275, 6.871, 7.250]),
            "withF1Constraint" : True 
        }
        ffpar.update({f"a_g_{iord}"   : 0.0 for iord in range(self.n["g"])})
        ffpar.update({f"a_f_{iord}"   : 0.0 for iord in range(self.n["f"])})
        ffpar.update({f"a_F1_{iord}"  : 0.0 for iord in range(self.n["F1"])})
        ffpar.update({f"a_F2_{iord}"  : 0.0 for iord in range(self.n["F2"])})
        ffpar.update({f"a_T1_{iord}"  : 0.0 for iord in range(self.n["T1"])})
        ffpar.update({f"a_T2_{iord}"  : 0.0 for iord in range(self.n["T2"])})
        ffpar.update({f"a_T23_{iord}" : 0.0 for iord in range(self.n["T23"])})
        return ffpar

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = np.sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr) # Beware arralike input here, will have problems

    def outer_fcns(self, q2: float) -> dict[str, float]:
        """Computation of the outer functions"""
        Mb = self.get_param(f"m_{self.B}")
        Mc = self.get_param(f"m_{self.V}")
        Mb2 = Mb*Mb
        Mb3 = Mb2*Mb
        rC = Mc/Mb
        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1.)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))

        nc    = self.get_userparam("nc")
        chif  = self.get_userparam("chif")
        chiF1 = self.get_userparam("chiF1")
        chiF2 = self.get_userparam("chiF2")
        chig  = self.get_userparam("chig")

        phif  = (4./Mb2)*np.sqrt(nc/(3*np.pi*chif))      * outer_function_vec(z, rC, 1, 3, 4)
        phiF1 = (4./Mb3)*np.sqrt(nc/(6*np.pi*chiF1))     * outer_function_vec(z, rC, 1, 5, 5)
        phiF2 = 8*np.sqrt(2)*np.sqrt(nc/(np.pi*chiF2))   * outer_function_vec(z, rC, 2, -1, 4)
        phig  = 16*np.sqrt(2)*np.sqrt(nc/(3*np.pi*chig)) * outer_function_vec(z, rC, 2, -1, 4)
        phif_0  = (4./Mb2)*np.sqrt(nc/(3*np.pi*chif))    * outer_function_vec(0.0, rC, 1, 3, 4)
        phiF1_0 = (4./Mb3)*np.sqrt(nc/(6*np.pi*chiF1))   * outer_function_vec(0.0, rC, 1, 5, 5)

        tp = (Mb+Mc)*(Mb+Mc)
        tm = (Mb-Mc)*(Mb-Mc)
        t0 = tm

        chiT1  = self.get_userparam("chiT1")
        chiT2  = self.get_userparam("chiT2")
        chiT23 = self.get_userparam("chiT23")
        phiT1  = np.sqrt(nc/(24*np.pi*chiT1))*outer_function_ten(q2, tp, tm, t0, 3, 3, 2)
        phiT2  = np.sqrt(nc*tp*tm/(24*np.pi*chiT2))*outer_function_ten(q2, tp, tm, t0, 1, 1, 2)
        phiT23 = np.sqrt(nc*((tp - tm)**2)/(48*tp*np.pi*chiT23))*outer_function_ten(q2, tp, tm, t0, 1, 1, 1)

        phi_fcn = {
            "Phig"    : phig   ,
            "Phif"    : phif   ,
            "PhiF1"   : phiF1  ,
            "PhiF2"   : phiF2  ,
            "Phif_0"  : phif_0 ,
            "PhiF1_0" : phiF1_0,
            "PhiT1"   : phiT1  ,
            "PhiT2"   : phiT2  ,
            "PhiT23"  : phiT23 ,
        }
        
        return phi_fcn

    def get_coef_arr(self, ffstr: str) -> list[float]:
        return np.array([self.get_userparam(f"a_{ffstr}_{iord}") for iord in range(self.n[ffstr])])

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Calculates BGL form factors following https://arxiv.org/pdf/2606.23410

        The choice is made for t0 = t- such that w=1 means z=0.
        
        Note there are kinematic constraints on tensor FFs in https://arxiv.org/pdf/2606.23410
        that are not applied here

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        Mb = self.get_param(f"m_{self.B}")
        Mc = self.get_param(f"m_{self.V}")
        
        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1.)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))
        zpow = np.array([z**ik for ik in range(int(self.n["nmax"]))])

        ag   = self.get_coef_arr("g")
        af   = self.get_coef_arr("f")
        aF1  = self.get_coef_arr("F1")
        aF2  = self.get_coef_arr("F2")
        aT1  = self.get_coef_arr("T1")
        aT2  = self.get_coef_arr("T2")
        aT23 = self.get_coef_arr("T23")
        
        # Blaschke factors
        Pf   = self.blaschke(self.get_userparam("BcStates1p"), z, Mb, Mc)
        PF1  = Pf
        Pg   = self.blaschke(self.get_userparam("BcStates1m"), z, Mb, Mc)
        PF2  = self.blaschke(self.get_userparam("BcStates0m"), z, Mb, Mc)
        PT1  = Pg
        PT2  = Pf
        PT23 = Pf
        
        # Outer functions
        outer_funcs = self.outer_fcns(q2)
        phig    = outer_funcs["Phig"]
        phif    = outer_funcs["Phif"]
        phiF1   = outer_funcs["PhiF1"]
        phif_0  = outer_funcs["Phif_0"]
        phiF1_0 = outer_funcs["PhiF1_0"]
        phiF2   = outer_funcs["PhiF2"]
        phiT1   = outer_funcs["PhiT1"]
        phiT2   = outer_funcs["PhiT2"]
        phiT23  = outer_funcs["PhiT23"]

        g = np.dot(ag, zpow[:len(ag)])/(Pg*phig)
        f = np.dot(af, zpow[:len(af)])/(Pf*phif)
        # https://arxiv.org/pdf/2606.23410, Eq. A.4, kinematic endpoint F1(1) = (Mb - Mc)f(1). 
        # - In C.23 it seems they use for zeorth order a_f = a_F2 * phi_f(z=0)/phi_F1(z=0) which I would assume is missing the (Mb - Mc)?
        F1 = 0
        if self.get_userparam("withF1Constraint"):
            F1 = (af[0]*(Mb-Mc)*phiF1_0/phif_0 + np.dot(aF1[1:], zpow[1:len(aF1)]))/(PF1*phiF1)
        else:
            F1 = np.dot(aF1, zpow[:len(aF1)])/(PF1*phiF1)
        # P1 = np.dot(aP1, zpow[:len(aP1)])*(np.sqrt(Mc/Mb)/((1+Mc/Mb)*PP1*phiP1))
        F2  = np.dot(aF2 , zpow[:len(aF2)])/(PF2*phiF2)
        T1  = np.dot(aT1 , zpow[:len(aT1)])/(PT1*phiT1)
        T2  = np.dot(aT2 , zpow[:len(aT2)])/(PT2*phiT2)
        T23 = np.dot(aT23, zpow[:len(aT23)])/(PT23*phiT23)

        # Relations taken from https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997-impl.hh lines 196 onwards
        ff = {
            "V"   : (Mb + Mc)*0.5*g,
            "A0"  : 0.5*F2         ,
            "A1"  : f/(Mb + Mc)    ,
            "A12" : F1/(8*Mb*Mc)   ,
            "T1"  : T1             ,
            "T2"  : T2             ,
            "T23" : T23
        }

        return ff




class FFBToV_BGL(FFBToV_BGLGeneric):
    _name = "FFBToV@BGL"
    def __init__(self, B: str, V: str):
        """BGL 2nd order in all FFs. Defualts are set from Type-A Baysian
        combined fit in https://arxiv.org/pdf/2606.23410 (Table 1 + Table 5 for tensor FFs)
        """
        super().__init__(B, V, 3, 3, 3, 3, 3, 3, 3)

    def define_userparams(self) -> dict[str, Any]:
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f_0"   : 0.01221 ,
            "a_f_1"   : 0.0135  , 
            "a_f_2"   : -0.16   ,
            "a_g_0"   : 0.03101 ,
            "a_g_1"   : -0.061  ,
            "a_g_2"   : -0.12   ,
            "a_F1_0"  : 0.002047, 
            "a_F1_1"  : -0.0038 , 
            "a_F1_2"  : -0.018  ,
            "a_F2_0"  : 0.05046 ,
            "a_F2_1"  : -0.210  ,
            "a_F2_2"  : 0.00    ,
            "a_T1_0"  : 0.01355 , 
            "a_T1_1"  : -0.034  , 
            "a_T1_2"  : -0.13   ,
            "a_T2_0"  : 0.003538, 
            "a_T2_1"  : 0.0004  , 
            "a_T2_2"  : -0.09   ,
            "a_T23_0" : 0.01055 , 
            "a_T23_1" : 0.000   , 
            "a_T23_2" : 0.13
        })


