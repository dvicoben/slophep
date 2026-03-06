import numpy as np

from slophep.Predictions.Observables import ObservableBase
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD2st

import flavio
from flavio.physics.bdecays.wilsoncoefficients import get_wceff_fccc_std
from flavio.physics import ckm

class BToD2stEllNuPrediction(ObservableBase):
    def __init__(self,
                 B: str,
                 M: str,
                 qiqj: str,
                 lep: str,
                 nu: str,
                 FF: FormFactorBToD2st,
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
    
    def _rate_prefactor(self, q2: float) -> float:
        # Prefactor Gamma_0 from Eq. (29) in https://arxiv.org/pdf/1711.03110
        GF = self.par['GF']
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self.B]
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.par)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.par)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)
        return p

    def dGdq2_SM(self, q2: float) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 32a
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.par["m_"+self.B]
        mC = self.par["m_"+self.M]
        r = mC/mB
        rhol = (self.par['m_'+self.lep]**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff(q2)
        kA1 = ff["kA1"]
        kA2 = ff["kA2"]
        kA3 = ff["kA3"]
        kV  = ff["kV"]
        wsqm1 = (w*w - 1)
        wmr = w-r
        wmrsq = wmr*wmr

        gamma = (2./3.)*gamma0*(r**3)*(wsqm1**1.5)*(q2hat - rhol)**2/(q2hat**3)*(
            kA1*kA1*(
                2*q2hat*(2*wmrsq + 3*q2hat)
                +rhol*(8*wmrsq - 3*q2hat)
            )
            + 2*wsqm1*(
                kA2*kA2*(
                    2*r*r*q2hat*wsqm1
                    +rhol*(3*q2hat + 4*r*r*wsqm1)
                )
                + kA3*kA3*(
                    2*q2hat*wsqm1
                    +rhol*(4*wmrsq - q2hat)
                )
                + 3*kV*kV*q2hat*(q2hat + 0.5*rhol)
                + 2*kA1*kA2*(
                    2*r*q2hat*wmr
                    +rhol*(3 - r*r - 2*r*w)
                )
                + 4*kA1*kA3*wmr*(q2hat + 2*rhol)
                + 2*kA2*kA3*(
                    2*r*q2hat*wsqm1
                    +rhol*(3*w*q2hat + 4*r*wsqm1)
                )
            )
        )
        return gamma


    def dGdq2(self, q2: float) -> float:
        # mb = running.get_mb(self.par, self.scale)
        # wc = get_wceff_fccc_std(self.wc_obj, self.par, self._qiqj, self.lep, self.nu, mb, self.scale, nf=5)
        # if self.lep != self.nu and all(C == 0 for C in wc.values()):
        #     # if all WCs vanish, so does the AC!
        #     return 0.0

        # SM only for now
        return self.dGdq2_SM(q2)

    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)
        
