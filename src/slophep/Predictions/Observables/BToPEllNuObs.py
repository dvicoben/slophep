from math import sqrt

from slophep.Predictions.FormFactorsBToP import FormFactorBToP
from slophep.Predictions.Observables import ObservableBase

import flavio
from flavio.physics.running import running
from flavio.physics.bdecays.wilsoncoefficients import get_wceff_fccc
from flavio.physics.bdecays import angular
from flavio.physics import ckm

class BToPEllNuPrediction(ObservableBase):
    def __init__(self, 
                 B: str,
                 P: str,
                 qiqj: str,
                 lep: str,
                 nu: str,
                 FF: FormFactorBToP,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        super().__init__(FF, ffargs, par, scale)
        self._B: str = B
        self._P: str = P
        self._qiqj: str = qiqj
        self._lep: str = lep
        self._nu: str = nu
        self.obslist = ['a', 'b', 'c']

    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def P(self) -> str:
        """The P meson"""
        return self._P
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (nu/e/tau)"""
        return self._nu
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
        mP = self.par['m_'+self._P]
        q2max = (mB-mP)**2
        return q2max
    
    def _prefactor(self, q2: float) -> float:
        """Return the prefactor including constants and CKM elements. Direct reimplementation
        of flavio equivalent in https://flav-io.github.io/apidoc/flavio/physics/bdecays/bplnu.m.html
        """
        GF = self.par['GF']
        ml = self.par['m_'+self.lep]
        # scale = self.scale
        # alphaem = running.get_alpha(self.par, scale)['alpha_e']
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.par)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.par)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0.0
        return 4.0*GF/sqrt(2)*Vij

    def get_angularcoeff(self, q2: float) -> dict:
        """Calculate angular coefficients, flavio is used for these calculations, method essentially follows
        flavio.physics.bdecays.bplnu._get_angularcoeff from https://flav-io.github.io/apidoc/flavio/physics/bdecays/bplnu.m.html

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents {a, b, c}
        """
        mb = running.get_mb(self.par, self.scale)
        wc = get_wceff_fccc(self.wc_obj, self.par, self._qiqj, self.lep, self.nu, mb, self.scale, nf=5)
        if self.lep != self.nu and all(C == 0 for C in wc.values()):
            return {'a': 0, 'b': 0, 'c': 0}  # if all WCs vanish, so does the AC!
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self.B]
        mP = self.par['m_'+self.P]
        mlight = running.get_mc(self.par, self.scale) if self._qiqj == "bc" else 0.0 # Set mlight = m_u = 0 for up quark "bu"
        N = self._prefactor(q2)
        ff = self.FF.get_ff(q2)
        h = angular.helicity_amps_p(q2, mB, mP, mb, mlight, ml, 0.0, ff, wc, N)
        J = angular.angularcoeffs_general_p(h, q2, mB, mP, mb, mlight, ml, 0)
        return J
    
    def dJ(self, q2: float) -> dict:
        return self.get_angularcoeff(q2)

    def J(self, q2: float) -> float:
        dJ = self.dJ(q2)
        dG = 2 * (dJ['a'] + dJ['c']/3.)
        # Can lead to Nan if dG==0 i.e. outside kinematic range, so:
        if dG <= 0.0:
            return {k : 0 for k in dJ}
        return {k : dJ[k]/dG for k in dJ}

    def dGdq2(self, q2: float) -> float:
        """Caclulate q2 distriution

        Parameters
        ----------
        q2 : float

        Returns
        -------
        float
            dGamma/dq2 (up to normalisation)
        """
        J = self.dJ(q2)
        return 2 * (J['a'] + J['c']/3.)
    
    def dBRdq2(self, q2: float) -> float:
        """Calculate differential BR, dBR/dq2
        
        Parameters
        ----------
        q2 : float

        Returns
        -------
        float
            dBR/dq2
        """
        dGdq2 = self.dGdq2(q2)
        BR = self.par[f"tau_{self._B}"]*dGdq2
        if self.P == 'pi0':
            # factor of 1/2 for neutral pi due to pi = (uubar-ddbar)/sqrt(2)
            return 0.5*BR
        return BR
    
    def afb(self, q2: float) -> float:
        return self.J(q2)["b"]