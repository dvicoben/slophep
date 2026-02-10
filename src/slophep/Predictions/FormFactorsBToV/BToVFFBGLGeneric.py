from math import sqrt
import numpy as np
from slophep.Predictions.FormFactorsBToV import FormFactorBToV

class BGLGeneric_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, 
                 order_g: int = 2, order_f: int = 2, order_F1: int = 2, order_F2: int = 1):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_BGLGeneric"
        self._ffpar = {}
        self._params = []
        self.params_setup(order_g, order_f, order_F1, order_F2)

        nmax = int(np.max([order_g+1, order_f+1, order_F1+1, order_F2+1])+1)

        internalparams = {
            "N_zexp_g"      : order_g+1,
            "N_zexp_f"      : order_f+1,
            "N_zexp_F1"     : order_F1+1,
            "N_zexp_F2"     : order_F2+1,
            "nmax"          : nmax,
            "Vcb"           : 41.5e-3,                       
            "chip"          : 5.131e-4,
            "chim"          : 3.894e-4,
            "chimL"         : 1.9421e-2,
            "BcStatesf"     : np.array([6.739, 6.750, 7.145, 7.150]),
            "BcStatesg"     : np.array([6.329, 6.920, 7.020]),
            "BcStatesP1"    : np.array([6.275, 6.842, 7.250]),
            "nc"            : 2.6,
            "etaEW"         : 1.0066,
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} Tensor FFs are 0, SM only parameterisation")

    def params_setup(self, order_g: int, order_f: int, order_F1: int, order_F2: int) -> None:
        self._params.extend([f"a{iord}" for iord in range(order_g + 1)])
        self._params.extend([f"b{iord}" for iord in range(order_f + 1)])
        self._params.extend([f"c{iord}" for iord in range(1, order_F1 + 1)])
        self._params.extend([f"d{iord}" for iord in range(order_F2 + 1)])
        
        self._ffpar = {k : 0.0 for k in self._params}

    def blaschke(self, BcStates: list, z: float, Mb: float, Mc: float) -> float:
        """Calculate Blaschke factor P(t) from B_c-type resonances"""
        Mb2 = Mb*Mb
        tp = (Mb+Mc)*(Mb+Mc)/Mb2
        tm = (Mb-Mc)*(Mb-Mc)/Mb2
        sqtptm = sqrt(tp - tm)
        sqtpBc = np.sqrt(tp-(BcStates/Mb)**2)
        parr = ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))))
        return np.prod(parr)

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

        phig = sqrt(256.*nc/(3*np.pi*chip))*((rC2*pow(1+z,2)*pow(1-z,-0.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phif = (1./Mb2)*sqrt(16.*nc/(3*np.pi*chim))*((rC*(1+z)*pow(1-z,1.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phiF1 = (1./Mb3)*sqrt(8.*nc/(3*np.pi*chim))*((rC*(1+z)*pow(1-z,2.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),5))
        phif_0 = 4.*rC*sqrt(nc/chim)/(Mb2*sqrt(3*np.pi)*pow(1+2*sqrC+rC,4))
        phiF1_0 = 2.*sqrt(2/(3*np.pi))*rC*sqrt(nc/chim)/(Mb3*pow(1+2*sqrC+rC,5))
        phiP1 = sqrt(nc/(np.pi*chimL))*((8.*sqrt(2)*rC2*pow(1+z,2)*pow(1-z,-0.5)))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4)

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
        z = (sqrt(w+1) - sqrt(2))/(sqrt(w+1) + sqrt(2))
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
        # return {k : ff[k]/(self.internalparams["etaEW"]*self.internalparams["Vcb"]) for k in ff} # This is what hammer seems to do
