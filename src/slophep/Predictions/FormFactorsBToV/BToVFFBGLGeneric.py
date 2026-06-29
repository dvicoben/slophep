import numpy as np
from slophep.Predictions.FormFactorsBToV import FormFactorBToV

class BGLGeneric_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, 
                 N_g: int = 3, N_f: int = 3, N_F1: int = 3, N_F2: int = 2):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_BGLGeneric"
        self._ffpar = {}
        self._params = []
        self.params_setup(N_g, N_f, N_F1, N_F2)

        nmax = int(np.max([N_g, N_f, N_F1, N_F2]))

        internalparams = {
            "N_zexp_g"      : N_g      ,
            "N_zexp_f"      : N_f      ,
            "N_zexp_F1"     : N_F1     ,
            "N_zexp_F2"     : N_F2     ,
            "nmax"          : nmax     ,
            "Vcb"           : 41.5e-3  ,                       
            "chip"          : 5.131e-4 ,
            "chim"          : 3.894e-4 ,
            "chimL"         : 1.9421e-2,
            "nc"            : 2.6      ,
            "etaEW"         : 1.0066   ,
            "BcStatesf"     : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStatesg"     : np.array([6.329, 6.920, 7.020]),
            "BcStatesP1"    : np.array([6.275, 6.842, 7.250]),
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} Tensor FFs are 0, SM only parameterisation")

    def params_setup(self, N_g: int, N_f: int, N_F1: int, N_F2: int) -> None:
        self._params.extend([f"a{iord}" for iord in range(N_g)])
        self._params.extend([f"b{iord}" for iord in range(N_f)])
        self._params.extend([f"c{iord}" for iord in range(1, N_F1)])
        self._params.extend([f"d{iord}" for iord in range(N_F2)])
        
        self._ffpar = {k : 0.0 for k in self._params}

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = np.sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr) # Beware with arraylike inputs here!

    def outer_fcns(self, z: float) -> dict:
        """Computation of the outer functions"""
        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]
        Mb2 = Mb*Mb
        Mb3 = Mb2*Mb
        rC = Mc/Mb
        rC2 = rC*rC
        sqrC = np.sqrt(rC)

        chim = self.internalparams["chim"]    # GeV^-2
        chip = self.internalparams["chip"]    # GeV^-2
        chimL = self.internalparams["chimL"]
        nc = self.internalparams["nc"]

        phig = np.sqrt(256.*nc/(3*np.pi*chip))*((rC2*np.power(1+z,2)*np.power(1-z,-0.5))/np.power(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phif = (1./Mb2)*np.sqrt(16.*nc/(3*np.pi*chim))*((rC*(1+z)*np.power(1-z,1.5))/np.power(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phiF1 = (1./Mb3)*np.sqrt(8.*nc/(3*np.pi*chim))*((rC*(1+z)*np.power(1-z,2.5))/np.power(((1+rC)*(1-z)+2*sqrC*(1+z)),5))
        phif_0 = 4.*rC*np.sqrt(nc/chim)/(Mb2*np.sqrt(3*np.pi)*np.power(1+2*sqrC+rC,4))
        phiF1_0 = 2.*np.sqrt(2/(3*np.pi))*rC*np.sqrt(nc/chim)/(Mb3*np.power(1+2*sqrC+rC,5))
        phiP1 = np.sqrt(nc/(np.pi*chimL))*((8.*np.sqrt(2)*rC2*np.power(1+z,2)*np.power(1-z,-0.5)))/np.power(((1+rC)*(1-z)+2*sqrC*(1+z)),4)

        phi_fcn = {
            "Phig" : phig,
            "Phif" : phif,
            "PhiF1" : phiF1,
            "PhiF2" : phiP1,
            "Phif_0" : phif_0,
            "PhiF1_0" : phiF1_0
        }
        
        return phi_fcn

    def get_ff(self, q2: float) -> dict:
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags
        
        Note that there is an additional 1./(etaEW*Vcb) factor applied to FFs in Hammer that is not
        present here. For that exact correspondence use BGL_Hammer FFs for the decay mode.

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
        
        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))
        zpow = np.array([z**ik for ik in range(int(self.internalparams["nmax"]))])

        ag  = np.array([self.ffpar[f"a{iord}"] for iord in range(self.internalparams["N_zexp_g"])])
        af  = np.array([self.ffpar[f"b{iord}"] for iord in range(self.internalparams["N_zexp_f"])])
        aF1 = np.array([self.ffpar[f"c{iord}"] for iord in range(1, self.internalparams["N_zexp_F1"])])
        aP1 = np.array([self.ffpar[f"d{iord}"] for iord in range(self.internalparams["N_zexp_F2"])])
        
        # Blaschke factors
        Pf = self.blaschke(self.internalparams["BcStatesf"], z, Mb, Mc)
        PF1 = Pf
        Pg = self.blaschke(self.internalparams["BcStatesg"], z, Mb, Mc)
        PP1 = self.blaschke(self.internalparams["BcStatesP1"], z, Mb, Mc)
        
        # Outer functions
        outer_funcs = self.outer_fcns(z)
        phig = outer_funcs["Phig"]
        phif = outer_funcs["Phif"]
        phiF1 = outer_funcs["PhiF1"]
        phif_0 = outer_funcs["Phif_0"]
        phiF1_0 = outer_funcs["PhiF1_0"]
        phiP1 = outer_funcs["PhiF2"]

        g = np.dot(ag, zpow[:len(ag)])/(Pg*phig)
        f = np.dot(af, zpow[:len(af)])/(Pf*phif)
        F1 = (af[0]*(Mb-Mc)*phiF1_0/phif_0 + np.dot(aF1, zpow[1:(1+len(aF1))]))/(PF1*phiF1)
        # P1 = np.dot(aP1, zpow[:len(aP1)])*(sqrt(Mc/Mb)/((1+Mc/Mb)*PP1*phiP1))
        F2 = np.dot(aP1, zpow[:len(aP1)])/(PP1*phiP1)

        # Relations taken from https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997-impl.hh lines 196 onwards
        ff = {
            "V"   : (Mb + Mc) * g / 2.0,
            "A0"  : 0.5*F2,
            "A1"  : f/(Mb + Mc),
            "A12" : F1/(8*Mb*Mc),
            "T1"  : 0.0,
            "T2"  : 0.0,
            "T3"  : 0.0,
            "T23" : 0.0
        }

        return ff



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



class BGLGeneric_BToV_wT(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, 
                 N_g: int = 3, N_f: int = 3, N_F1: int = 3, N_F2: int = 2,
                 N_T1: int = 0, N_T2: int = 0, N_T23: int = 0):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_BGLGeneric"
        self._ffpar = {}
        self._params = []
        self.params_setup(N_g, N_f, N_F1, N_F2, N_T1, N_T2, N_T23)

        nmax = int(np.max([N_g, N_f, N_F1, N_F2, N_T1, N_T2, N_T23]))

        internalparams = {
            "N_zexp_g"      : N_g      ,
            "N_zexp_f"      : N_f      ,
            "N_zexp_F1"     : N_F1     ,
            "N_zexp_F2"     : N_F2     ,
            "N_zexp_T1"     : N_T1     ,
            "N_zexp_T2"     : N_T2     ,
            "N_zexp_T23"    : N_T23    ,
            "nmax"          : nmax     ,
            "Vcb"           : 41.5e-3  ,                       
            "chig"          : 5.131e-4 ,
            "chif"          : 3.894e-4 ,
            "chiF1"         : 3.894e-4 ,
            "chiF2"         : 1.9421e-2,
            "chiT1"         : 4.98e-4  ,
            "chiT2"         : 2.77e-4  ,
            "chiT23"        : 2.77e-4  ,
            "nc"            : 2.6      ,
            "etaEW"         : 1.0066   ,
            "BcStates1p"    : np.array([6.739, 6.750, 7.145, 7.150]), #A1
            "BcStates1m"    : np.array([6.329, 6.910, 7.020]),        #V1
            "BcStates0m"    : np.array([6.275, 6.871, 7.250]),        #A0
        }
        self._internalparams.update(internalparams)

    def params_setup(self, N_g: int, N_f: int, N_F1: int, N_F2: int, N_T1: int, N_T2: int, N_T23: int) -> None:
        self._params.extend([f"g_{iord}"   for iord in range(N_g)    ])
        self._params.extend([f"f_{iord}"   for iord in range(N_f)    ])
        self._params.extend([f"F1_{iord}"  for iord in range(1, N_F1)])
        self._params.extend([f"F2_{iord}"  for iord in range(N_F2)   ])
        self._params.extend([f"T1_{iord}"  for iord in range(N_T1)   ])
        self._params.extend([f"T2_{iord}"  for iord in range(N_T2)   ])
        self._params.extend([f"T13_{iord}" for iord in range(N_T23)  ])
        
        self._ffpar = {k : 0.0 for k in self._params}

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = np.sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr) # Beware arralike input here, will have problems

    def outer_fcns(self, q2: float) -> dict:
        """Computation of the outer functions"""
        Mb = self.internalparams["Mb"]
        Mc = self.internalparams["Mc"]
        Mb2 = Mb*Mb
        Mb3 = Mb2*Mb
        rC = Mc/Mb
        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1.)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))

        nc    = self.internalparams["nc"]
        chif  = self.internalparams["chif"]
        chiF1 = self.internalparams["chiF1"]
        chiF2 = self.internalparams["chiF2"]
        chig  = self.internalparams["chig"]

        phif  = (4./Mb2)*np.sqrt(nc/(3*np.pi*chif))      * outer_function_vec(z, rC, 1, 3, 4)
        phiF1 = (4./Mb3)*np.sqrt(nc/(6*np.pi*chiF1))     * outer_function_vec(z, rC, 1, 5, 5)
        phiF2 = 8*np.sqrt(2)*np.sqrt(nc/(np.pi*chiF2))   * outer_function_vec(z, rC, 2, -1, 4)
        phig  = 16*np.sqrt(2)*np.sqrt(nc/(3*np.pi*chig)) * outer_function_vec(z, rC, 2, -1, 4)
        phif_0  = (4./Mb2)*np.sqrt(nc/(3*np.pi*chif))    * outer_function_vec(0.0, rC, 1, 3, 4)
        phiF1_0 = (4./Mb3)*np.sqrt(nc/(6*np.pi*chiF1))   * outer_function_vec(0.0, rC, 1, 5, 5)

        tp = (Mb+Mc)*(Mb+Mc)
        tm = (Mb-Mc)*(Mb-Mc)
        t0 = tm

        chiT1  = self.internalparams["chiT1"]
        chiT2  = self.internalparams["chiT2"]
        chiT23 = self.internalparams["chiT23"]
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

    def get_ff(self, q2: float) -> dict:
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags
        
        Note that there is an additional 1./(etaEW*Vcb) factor applied to FFs in Hammer that is not
        present here. For that exact correspondence use BGL_Hammer FFs for the decay mode.

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
        
        w = max((Mb**2 + Mc**2 - q2) / (2 * Mb * Mc), 1.)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))
        zpow = np.array([z**ik for ik in range(int(self.internalparams["nmax"]))])

        ag   = np.array([self.ffpar[f"g_{iord}"]   for iord in range(self.internalparams["N_zexp_g"])])
        af   = np.array([self.ffpar[f"f_{iord}"]   for iord in range(self.internalparams["N_zexp_f"])])
        aF1  = np.array([self.ffpar[f"F1_{iord}"]  for iord in range(1, self.internalparams["N_zexp_F1"])])
        aF2  = np.array([self.ffpar[f"F2_{iord}"]  for iord in range(self.internalparams["N_zexp_F2"])])
        aT1  = np.array([self.ffpar[f"T1_{iord}"]  for iord in range(self.internalparams["N_zexp_T1"])])
        aT2  = np.array([self.ffpar[f"T2_{iord}"]  for iord in range(self.internalparams["N_zexp_T2"])])
        aT23 = np.array([self.ffpar[f"T23_{iord}"] for iord in range(self.internalparams["N_zexp_T23"])])
        
        # Blaschke factors
        Pf   = self.blaschke(self.internalparams["BcStates1p"], z, Mb, Mc)
        PF1  = Pf
        Pg   = self.blaschke(self.internalparams["BcStates1m"], z, Mb, Mc)
        PF2  = self.blaschke(self.internalparams["BcStates0m"], z, Mb, Mc)
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
        F1 = (af[0]*(Mb-Mc)*phiF1_0/phif_0 + np.dot(aF1, zpow[1:(1+len(aF1))]))/(PF1*phiF1)
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
