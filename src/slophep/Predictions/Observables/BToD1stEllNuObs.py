import numpy as np

from slophep.Predictions.Observables import ObservableBase
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD1st
from slophep.Predictions.Math.wilson_coefs import get_wceff_fccc

import flavio
from flavio.physics import ckm

class BToD1stEllNuPrediction(ObservableBase):
    def __init__(self,
                 B: str,
                 M: str,
                 qiqj: str,
                 lep: str,
                 nu: str,
                 FF: FormFactorBToD1st,
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
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)/(192*(np.pi**3))
        return p

    def dGdq2M_SM(self, q2: float, mC: float = None) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 31a
        # And https://scoap3-prod-backend.s3.cern.ch/media/files/79991/10.1088/1674-1137/ace821.pdf eq. 109
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.par["m_"+self.B]
        mC = self.par["m_"+self.M] if mC is None else mC
        r = mC/mB
        rhol = (self.par['m_'+self.lep]**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        # ff = self.FF.get_ff(q2) if mC is None else self.FF.get_ff_mmeson(q2, mC)
        ff = self.FF.get_ff_mmeson(q2, mC)
        fV1 = ff["gV1"]
        fV2 = ff["gV2"]
        fV3 = ff["gV3"]
        fa  = ff["gA"]
        wsqm1 = (w**2 - 1)

        gamma = gamma0*(1/(mB*mC))*(r**3)*np.sqrt(wsqm1)*(q2hat - rhol)**2/(q2hat**3)*(
            fV1**2 * (
                2*q2hat*((w - r)**2 + 2*q2hat)
                + rhol*(4*((w - r)**2) - q2hat)
            )
            + wsqm1 * (
                fV2**2 * (
                    2*(r**2)*q2hat*wsqm1
                    + rhol*(3*q2hat + 4*(r**2)*wsqm1)
                )
                + fV3**2 * (
                    2*q2hat*wsqm1
                    + rhol*(4*((w - r)**2) - q2hat)
                )
                + 2*(fa**2)*q2hat*(2*q2hat + rhol)
                + 2*fV1*fV2 * (
                    2*r*q2hat*(w-r)
                    + rhol*(3 - r**2 - 2*r*w)
                )
                + 4*fV1*fV3*(w-r)*(q2hat + 2*rhol)
                + 2*fV2*fV3*(
                    2*r*q2hat*wsqm1
                    + rhol*(3*w*q2hat + 4*r*wsqm1)
                )
            )
        )
        return gamma

    def dGdq2_SM(self, q2: float) -> float:
        return self.dGdq2M_SM(q2, None)

    def dGdq2M(self, q2: float, mC: float) -> float:
        gammaSM = self.dGdq2M_SM(q2, mC)
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.

        wc = get_wceff_fccc(self.wc_obj, self.par, self._qiqj, self.lep, self.nu, self.scale, nf=5, withSM = False)

        if self.lep != self.nu and all(C == 0 for C in wc.values()):
            return 0.0
        # Early return if no NP
        if all(C == 0 for C in wc.values()):
            return gammaSM
        
        Cvl = wc["VL"]
        Cvr = wc["VR"]
        Csl = wc["SL"]
        Csr = wc["SR"]
        Ct  = wc["T"]

        mB = self.par["m_"+self.B]
        mC = self.par["m_"+self.M] if mC is None else mC
        r = mC/mB
        rhol = (self.par['m_'+self.lep]**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC)
        fV1 = ff["gV1"]
        fV2 = ff["gV2"]
        fV3 = ff["gV3"]
        fT1 = ff["gT1"]
        fT2 = ff["gT2"]
        fT3 = ff["gT3"]
        fa  = ff["gA"]
        fs  = ff["gS"]
        wsqm1 = (w**2 - 1)

        # This should follow https://arxiv.org/pdf/1711.03110 Eq. 31b, needs to be checked, high risk of typos
        nscale = gamma0*(1/(mB*mC))*(r**3)*np.sqrt(wsqm1)*(q2hat - rhol)**2/(q2hat**3)
        gamma = gammaSM*(np.abs(1. + Cvl - Cvr)**2) + nscale*(
            3*fs*fs*wsqm1*(q2hat**2)*np.abs(Csl + Csr)**2
            - 6*fs*wsqm1*np.sqrt(rhol)*q2hat*(np.real(Csr-Csl) 
                 + np.real((Csr-Csl)*np.conjugate(Cvl)) 
                 + np.real((Csr-Csl)*np.conjugate(Cvr))
            )*(
                fV1 + fV2*(1 - r*w) + fV3*(w-r)
            )
            + 16*(np.abs(Ct)**2)*(q2hat + 2*rhol)*(
                fT1*fT1*(
                    q2hat*(2+w*w)
                    + 4*r*r*wsqm1
                )
                + fT2*fT2*(w * (w-r)**2 - q2hat)
                + fT3*fT3*q2hat*wsqm1*wsqm1
                + 2*fT1*fT2*(3*w*q2hat + 4*r*wsqm1) 
                - 2*fT3*(fT1*w + fT2)*q2hat*wsqm1
            )
            - 24*np.sqrt(rhol)*q2hat*(
                (np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) - np.real(Ct*np.conjugate(Cvr)))*(
                    2*fa*(fT1*r + fT2)*wsqm1
                )
                - (np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) + np.real(Ct*np.conjugate(Cvr)))*(
                    2*fT1*fV1*(1 - r*w) + (w*fT1 + 3*fT2 - fT3*wsqm1)*fV1*(w -r)
                    + (w*fT1 + fT2 - fT3*wsqm1)*(fV2*r + fV3)*wsqm1
                )
            )
            + 4*(np.real(Cvr) + np.real(Cvr*np.conjugate(Cvl)))*(
                3*fV1*fV1*q2hat(2*q2hat + rhol) 
                + 2*fV1*(fV1 + 2*fV3*(w - r))*wsqm1*(q2hat + 2*rhol)
                + wsqm1*(
                    fV2*fV2*(2*r*r*q2hat*wsqm1 + rhol*(3*q2hat + 4*r*r*wsqm1)) 
                    + fV3*fV3*(2*q2hat*wsqm1 + rhol(4*(w - r)**2 - q2hat))
                    + 2*fV1*fV2*(2*r*q2hat*(w - r) + rhol*(3 - r*r - 2*r*w)) 
                    + 2*fV2*fV3*(2*r*q2hat*wsqm1 + rhol*(3*w*q2hat + 4*r*wsqm1))
                )
            )
        )
        return gamma
    
    def dGdq2(self, q2: float) -> float:
        return self.dGdq2M(q2, None)

    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)
        
