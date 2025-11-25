import slophep.Predictions.Math.BaryonicMathTools as bmt
from slophep.Predictions.Observables import ObservableBase
from slophep.Predictions.FormFactorsBaryonic import FormFactorOneHalfpToOneHalfp

from flavio.physics.running import running
from flavio.physics.bdecays.wilsoncoefficients import get_wceff_fccc_std


class LbToLcEllNuPredictionRaw(ObservableBase):
    def __init__(self,
                 B: str,
                 M: str,
                 qiqj: str,
                 lep: str,
                 nu: str,
                 FF: FormFactorOneHalfpToOneHalfp,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8
                 ):
        super().__init__(FF, ffargs, par, scale)
        self._B: str = B
        self._M: str = M
        self._qiqj: str = qiqj
        self._lep: str = lep
        self._nu: str = nu
        self.obslist = ["1ss", "1cc", "1c", "2ss", "2cc", "2c", "3sc", "3s", "4sc", "4s"]
        self._alphaL = -0.768

    @property
    def B(self) -> str:
        """The B baryon"""
        return self._B
    @property
    def M(self) -> str:
        """The final baryon"""
        return self._M
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (nu/e/tau)"""
        return self._nu
    @property
    def alphaL(self) -> float:
        """Lambda decay constant"""
        return self._alphaL
    @property
    def q2min(self) -> float:
        """Minimum physical q2"""
        ml = self.par['m_'+self.lep]
        q2min = ml**2
        return q2min
    @property
    def q2max(self) -> float:
        """Maximum physical q2"""
        mB = self.par['m_'+self._B]
        mM = self.par['m_'+self._M]
        q2max = (mB-mM)**2
        return q2max
    
    def set_alphaL(self, alpha: float):
        """Set Lambda decay constant"""
        self._alphaL = alpha

    def _prefactor(self, q2: float) -> float:
        mB = self.par['m_'+self._B]
        mM = self.par['m_'+self._M]
        ml = self.par['m_'+self.lep]
        return bmt.LbToLcEllNu_norm_EOS(q2, mB, mM, ml, self.par, self._qiqj)
    
    def get_angularcoeff(self, q2: float) -> dict:
        """Calculate angular coefficients, as per https://arxiv.org/abs/1907.12554,
        reproducing the EOS implementation https://github.com/eos/eos/blob/v1.0.13/eos/b-decays/lambdab-to-lambdac-l-nu.cc

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents K_i
        """
        mb = running.get_mb(self.par, self.scale)
        wc = get_wceff_fccc_std(self.wc_obj, self.par, self._qiqj, self.lep, self.nu, mb, self.scale, nf=5)
        if self.lep != self.nu and all(C == 0 for C in wc.values()):
            # if all WCs vanish, so does the AC!
            return {k: 0 for k in self.obslist}
        
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self.B]
        mM = self.par['m_'+self.M]
        mlight = running.get_mc(self.par, self.scale) if self._qiqj == "bc" else 0.0 # Set mlight = m_u = 0 for up quark "bu"
        N = self._prefactor(q2)
        ff = self.FF.get_ff(q2)
        A = bmt.LbToLcEllNu_amplitudes(q2, mB, mM, mb, mlight, ml, ff, wc, N)
        K = bmt.LbToLcEllNu_observables(q2, mB, mM, ml, A, self.alphaL)
        return K
    
    def dJ(self, q2: float) -> dict:
        """Alias for get\_angularcoeff

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents K_i
        """
        return self.get_angularcoeff(q2)

    def J(self, q2: float) -> dict:
        """Rate-normalised angular coefficients

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents K_i
        """
        K = self.dJ(q2)
        norm = self._dGdq2(K)
        return {iobs : ik/norm for iobs, ik in K.items()}

    def _dGdq2(self, K: dict) -> float:
        return 2.0*K["1ss"] + K["1cc"]

    def dGdq2(self, q2: float) -> float:
        """Decay rate over q2
        
        Parameters
        ----------
        q2 : float

        Returns
        -------
        float
            dGamma/dq2 (up to normalisation)
        """
        K = self.dJ(q2)
        return self._dGdq2(K)
    
    def dBRdq2(self, q2: float) -> float:
        """Differential branching ratio

        Parameters
        ----------
        q2 : float

        Returns
        -------
        float
            dBR/dq2
        """
        dG = self.dGdq2(q2)
        BR = dG*self.par[f"tau_{self.B}"]
        # tau = 1.471e-12
        # hbar = 6.58211928e-25
        # BR = dG*(tau/hbar)
        return BR
    
    def afb_lep(self, q2: float) -> float:
        """Leptonic Afb"""
        K = self.dJ(q2)
        dG = self._dGdq2(K)
        return 3.0 / 2.0 * K["1c"] / dG
    
    def afb_had(self, q2: float) -> float:
        """Hadronic Afb"""
        K = self.dJ(q2)
        dG = self._dGdq2(K)
        return 1.0 / 2.0 * (2.0 * K["2ss"] + K["2cc"]) / dG
    
    def afb_comb(self, q2: float) -> float:
        """Combined Afb"""
        K = self.dJ(q2)
        dG = self._dGdq2(K)
        return 3.0 / 4.0 * K["2c"] / dG
    
    def f0(self, q2: float) -> float:
        """F0 Observable"""
        K = self.dJ(q2)
        dG = self._dGdq2(K)
        return (2.0 * K["1ss"] - K["1cc"]) / dG