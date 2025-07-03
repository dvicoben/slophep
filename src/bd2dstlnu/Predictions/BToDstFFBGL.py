from math import sqrt
import numpy as np
from bd2dstlnu.Predictions.BToDstFFBase import FormFactor


class BGL(FormFactor):
    def __init__(self, par: dict = None, scale: float = None, *ffargs):
        super().__init__(par, scale)
        
        self._name = "BGL"
        self._ffpar = { # https://arxiv.org/pdf/1707.09509 (Third column in table V)
            "a0" : 0.0209,
            "a1" : 0.33,
            "a2" : 0.6,
            "b0" : 0.01218,
            "b1" : 0.046,
            "b2" : 0.48,
            "c1" : 0.0063,
            "c2" : 0.062,
            "d0" : 0.0595,
            "d1" : -0.218
        }
        self._params = ["a0", "a1", "a2", "b0", "b1", "b2", "c1", "c2", "d0", "d1"]
        self._process = 'B->D*'
        self._pd = {'B': 'B0', 'V': 'D*+', 'q': 'b->c'}

        self._internalparams = {
            "Mb"         : self.par['m_'+self._pd['B']],
            "Mc"         : self.par['m_'+self._pd['V']],
            "nmax"       : 4,
            "Vcb"        : 41.5e-3,                       
            "chim"       : 3.068e-4, # GeV^-2
            "chip"       : 5.280e-4, # GeV^-2
            "chimL"      : 2.466e-3,
            "nc"         : 2.6,
            "etaEW"      : 1.0066,
            "BcStatesf"  : np.array([6.730,6.736,7.135,7.142]),  # GeV
            "BcStatesg"  : np.array([6.337,6.899,7.012,7.280]),  # GeV
            "BcStatesP1" : np.array([6.275,6.842,7.250])         # GeV
        }

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
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags

        Could alternatively consider EOS implementation
        https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997.hh
        https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997-impl.hh

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

        ag  = np.array([self.ffpar["a0"], self.ffpar["a1"], self.ffpar["a2"]])
        af  = np.array([self.ffpar["b0"], self.ffpar["b1"], self.ffpar["b2"]])
        aF1 = np.array([self.ffpar["c1"], self.ffpar["c2"]])
        aP1 = np.array([self.ffpar["d0"], self.ffpar["d1"]])

        chim = self.internalparams["chim"]    # GeV^-2
        chip = self.internalparams["chip"]    # GeV^-2
        chimL = self.internalparams["chimL"]
        nc = self.internalparams["nc"]

        # Blaschke factors
        Pf = self.blaschke(self.internalparams["BcStatesf"], z, Mb, Mc)
        PF1 = Pf
        Pg = self.blaschke(self.internalparams["BcStatesg"], z, Mb, Mc)
        PP1 = self.blaschke(self.internalparams["BcStatesP1"], z, Mb, Mc)
        
        phig = sqrt(256.*nc/(3*np.pi*chip))*((rC2*pow(1+z,2)*pow(1-z,-0.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phif = (1./Mb2)*sqrt(16.*nc/(3*np.pi*chim))*((rC*(1+z)*pow(1-z,1.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4))
        phiF1 = (1./Mb3)*sqrt(8.*nc/(3*np.pi*chim))*((rC*(1+z)*pow(1-z,2.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),5))
        phif_0 = 4.*rC*sqrt(nc/chim)/(Mb2*sqrt(3*np.pi)*pow(1+2*sqrC+rC,4))
        phiF1_0 = 2.*sqrt(2/(3*np.pi))*rC*sqrt(nc/chim)/(Mb3*pow(1+2*sqrC+rC,5))
        phiP1 = sqrt(nc/(np.pi*chimL))*((8.*sqrt(2)*rC2*pow(1+z,2)*pow(1-z,-0.5)))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4)

        g = np.dot(ag, zpow[:len(ag)])/(Pg*phig)
        f = np.dot(af, zpow[:len(af)])/(Pf*phif)
        F1 = (af[0]*(Mb-Mc)*phiF1_0/phif_0 + np.dot(aF1, zpow[1:(1+len(aF1))]))/(PF1*phiF1)
        # P1 = np.dot(aP1, zpow[:len(aP1)])*(sqrC/((1+rC)*PP1*phiP1))
        F2 = np.dot(aP1, zpow[:len(aP1)])/(PP1*phiP1)

        # Relations taken from https://github.com/eos/eos/blob/master/eos/form-factors/parametric-bgl1997-impl.hh lines 196 onwards
        ff = {
            "V"  : (Mb + Mc) * g / 2.0,
            "A0" : 0.5*F2,
            "A1" : f/(Mb + Mc),
            "A12" : F1/(8*Mb*Mc),
            "T1" : 0.0,
            "T2" : 0.0,
            "T23" : 0.0
        }

        return ff
        # return {k : ff[k]/(self.internalparams["etaEW"]*self.internalparams["Vcb"]) for k in ff} # This is what hammer seems to do
