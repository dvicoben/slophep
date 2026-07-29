from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

from slophep.FormFactors.FormFactorsBToDstst import (
    FormFactorBToD0st, FormFactorBToD1, FormFactorBToD1st, FormFactorBToD2st
)
from slophep.Observables.ObservableImpl import Observable
from slophep.Tools.SamplingTools import fluctsettings, FluctType

from slophep.Math.wilson_coefs import get_wceff_fccc

import flavio
from flavio.physics import ckm



class BToD0stEllNuPrediction(Observable):
    _name = "ObsBToD0stEllNu"
    def __init__(self, 
                 B   : str,
                 M   : str,
                 lep : str,
                 nu  : str,
                 FF  : FormFactorBToD0st):
        super().__init__(FF)
        self._B: str = B
        self._M: str = M
        self._qiqj: str = "bc"
        self._lep: str = lep
        self._nu: str = nu

    @property
    def scale(self) -> float:
        return self.pm.scale
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def M(self) -> str:
        """The final-state meson"""
        return self._M
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (mu/e/tau)"""
        return self._nu
    @property
    def q2min(self) -> float:
        """Minimum physical q2"""
        ml = self.get_param('m_'+self.lep)
        q2min = ml**2
        return q2min
    @property
    def q2max(self) -> float:
        """Maximum physical q2"""
        mB = self.get_param('m_'+self.B)
        mM = self.get_param('m_'+self.M)
        q2max = (mB-mM)**2
        return q2max

    def _rate_prefactor(self, q2: float) -> float:
        # Prefactor Gamma_0 from Eq. (29) in https://arxiv.org/pdf/1711.03110
        GF = self.get_param('GF')
        ml = self.get_param('m_'+self.lep)
        mB = self.get_param('m_'+self.B)
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.pm)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.pm)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)/(192*(np.pi**3))
        return p

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M_SM(self, q2: float, mC: float = None) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 31a
        # And https://scoap3-prod-backend.s3.cern.ch/media/files/79991/10.1088/1674-1137/ace821.pdf eq. 109
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC, mB)
        gp = ff["g+"]
        gm = ff["g-"]
        wsqm1 = (w**2 - 1)

        gamma = 2*gamma0*(1/(mB*mC))*(r**3)*np.sqrt(wsqm1)*(q2hat - rhol)**2/(q2hat**3)*(
            gm*gm*(w-1)*(
                rhol*((1+r*r)*(2*w-1) + 2*r*(w-2))
                + ((1-r)**2)*(w+1)*q2hat
            )
            + gp*gp*(w+1)*(
                rhol*((1+r*r)*(2*w+1) - 2*r*(w+2))
                + ((1+r)**2)*(w-1)*q2hat
            )
            -2*gm*gp*(1-r*r)*wsqm1*(q2hat + 2*rhol)
        )
        return gamma
    
    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M(self, q2: float, mC: float) -> float:
        gammaSM = self.dGdq2M_SM(q2, mC)
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.

        wc = get_wceff_fccc(self.pm.get_wc(), self.pm, self._qiqj, self.lep, self.nu, self.scale, nf=5, withSM = False)

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

        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC)
        gps = ff["gP"]
        gp  = ff["g+"]
        gm  = ff["g-"]
        gt  = ff["gT"]
        wsqm1 = (w**2 - 1)

        # This should follow https://arxiv.org/pdf/1711.03110 Eq. 30b, needs to be checked, high risk of typos
        nscale = gamma0*(1/(mB*mC))*(r**3)*np.sqrt(wsqm1)*(q2hat - rhol)**2/(q2hat**2)
        gamma = gammaSM*(np.abs(1. + Cvl - Cvr)**2) + nscale*(
            3*q2hat*gps*gps*np.abs(Csr - Csl)**2 
            + 6*(np.real(Csr-Csl) 
                    + np.real((Csr-Csl)*np.conjugate(Cvl)) 
                    - np.real((Csr-Csl)*np.conjugate(Cvr))
                )*gps*np.sqrt(rhol)*(
                gm*(1+r)*(w-1) 
                - gp*(1-r)*(w+1)
            )
            + 8*gt*wsqm1*(
                2*np.abs(Ct)*gt*(q2hat + 2*rhol)
                + 3*np.sqrt(rhol)*(
                    np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) - np.real(Ct*np.conjugate(Cvr))
                )*(
                    gp*(1+r)
                    - gm*(1-r)
                )
            )
        )
        return gamma

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2_SM(self, q2: float) -> float:
        return self.dGdq2M_SM(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2(self, q2: float) -> float:
        return self.dGdq2M(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)
        



class BToD1EllNuPrediction(Observable):
    _name = "ObsBToD1EllNu"
    def __init__(self, 
                 B   : str,
                 M   : str,
                 lep : str,
                 nu  : str,
                 FF  : FormFactorBToD1):
        super().__init__(FF)
        self._B: str = B
        self._M: str = M
        self._qiqj: str = "bc"
        self._lep: str = lep
        self._nu: str = nu

    @property
    def scale(self) -> float:
        return self.pm.scale
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def M(self) -> str:
        """The final-state meson"""
        return self._M
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (mu/e/tau)"""
        return self._nu
    @property
    def q2min(self) -> float:
        """Minimum physical q2"""
        ml = self.get_param('m_'+self.lep)
        q2min = ml**2
        return q2min
    @property
    def q2max(self) -> float:
        """Maximum physical q2"""
        mB = self.get_param('m_'+self.B)
        mM = self.get_param('m_'+self.M)
        q2max = (mB-mM)**2
        return q2max

    def _rate_prefactor(self, q2: float) -> float:
        # Prefactor Gamma_0 from Eq. (29) in https://arxiv.org/pdf/1711.03110
        GF = self.get_param('GF')
        ml = self.get_param('m_'+self.lep)
        mB = self.get_param('m_'+self.B)
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.pm)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.pm)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)/(192*(np.pi**3))
        return p

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M_SM(self, q2: float, mC: float = None) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 31a
        # And https://scoap3-prod-backend.s3.cern.ch/media/files/79991/10.1088/1674-1137/ace821.pdf eq. 109
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC, mB)
        fV1 = ff["fV1"]
        fV2 = ff["fV2"]
        fV3 = ff["fV3"]
        fa  = ff["fA"]
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
    
    @fluctsettings(FluctType.NUMERIC)
    def dGdq2_SM(self, q2) -> float:
        return self.dGdq2M_SM(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M(self, q2: float, mC: float) -> float:
        gammaSM = self.dGdq2M_SM(q2, mC)
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.

        wc = get_wceff_fccc(self.pm.get_wc(), self.pm, self._qiqj, self.lep, self.nu, self.scale, nf=5, withSM = False)

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

        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC)
        fV1 = ff["fV1"]
        fV2 = ff["fV2"]
        fV3 = ff["fV3"]
        fT1 = ff["fT1"]
        fT2 = ff["fT2"]
        fT3 = ff["fT3"]
        fa  = ff["fA"]
        fs  = ff["fS"]
        wsqm1 = (w**2 - 1)

        # This should follow https://arxiv.org/pdf/1711.03110 Eq. 31b, needs to be checked, high risk of typos
        nscale = gamma0*(1/(mB*mC))*(r**3)*np.sqrt(wsqm1)*(q2hat - rhol)**2/(q2hat**3)
        gamma = gammaSM*(np.abs(1. + Cvl - Cvr)**2) + nscale*(
            3*fs*fs*wsqm1*(q2hat**2)*np.abs(Csl + Csr)**2
            + 6*fs*wsqm1*np.sqrt(rhol)*q2hat*(np.real(Csr-Csl) 
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

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2(self, q2: float) -> float:
        return self.dGdq2M(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)


        
class BToD1stEllNuPrediction(Observable):
    _name = "ObsBToD1stEllNu"
    def __init__(self, 
                 B   : str,
                 M   : str,
                 lep : str,
                 nu  : str,
                 FF  : FormFactorBToD1st):
        super().__init__(FF)
        self._B: str = B
        self._M: str = M
        self._qiqj: str = "bc"
        self._lep: str = lep
        self._nu: str = nu

    @property
    def scale(self) -> float:
        return self.pm.scale
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def M(self) -> str:
        """The final-state meson"""
        return self._M
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (mu/e/tau)"""
        return self._nu
    @property
    def q2min(self) -> float:
        """Minimum physical q2"""
        ml = self.get_param('m_'+self.lep)
        q2min = ml**2
        return q2min
    @property
    def q2max(self) -> float:
        """Maximum physical q2"""
        mB = self.get_param('m_'+self.B)
        mM = self.get_param('m_'+self.M)
        q2max = (mB-mM)**2
        return q2max

    def _rate_prefactor(self, q2: float) -> float:
        # Prefactor Gamma_0 from Eq. (29) in https://arxiv.org/pdf/1711.03110
        GF = self.get_param('GF')
        ml = self.get_param('m_'+self.lep)
        mB = self.get_param('m_'+self.B)
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.pm)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.pm)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)/(192*(np.pi**3))
        return p

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M_SM(self, q2: float, mC: float = None) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 31a
        # And https://scoap3-prod-backend.s3.cern.ch/media/files/79991/10.1088/1674-1137/ace821.pdf eq. 109
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
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

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2_SM(self, q2: float) -> float:
        return self.dGdq2M_SM(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M(self, q2: float, mC: float) -> float:
        gammaSM = self.dGdq2M_SM(q2, mC)
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.

        wc = get_wceff_fccc(self.pm.get_wc(), self.pm, self._qiqj, self.lep, self.nu, self.scale, nf=5, withSM = False)

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

        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC, mB)
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
    
    @fluctsettings(FluctType.NUMERIC)
    def dGdq2(self, q2: float) -> float:
        return self.dGdq2M(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)
        


        
class BToD2stEllNuPrediction(Observable):
    _name = "ObsBToD2stEllNu"
    def __init__(self, 
                 B   : str,
                 M   : str,
                 lep : str,
                 nu  : str,
                 FF  : FormFactorBToD2st):
        super().__init__(FF)
        self._B: str = B
        self._M: str = M
        self._qiqj: str = "bc"
        self._lep: str = lep
        self._nu: str = nu

    @property
    def scale(self) -> float:
        return self.pm.scale
    @property
    def B(self) -> str:
        """The B meson"""
        return self._B
    @property
    def M(self) -> str:
        """The final-state meson"""
        return self._M
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (mu/e/tau)"""
        return self._nu
    @property
    def q2min(self) -> float:
        """Minimum physical q2"""
        ml = self.get_param('m_'+self.lep)
        q2min = ml**2
        return q2min
    @property
    def q2max(self) -> float:
        """Maximum physical q2"""
        mB = self.get_param('m_'+self.B)
        mM = self.get_param('m_'+self.M)
        q2max = (mB-mM)**2
        return q2max

    def _rate_prefactor(self, q2: float) -> float:
        # Prefactor Gamma_0 from Eq. (29) in https://arxiv.org/pdf/1711.03110
        GF = self.get_param('GF')
        ml = self.get_param('m_'+self.lep)
        mB = self.get_param('m_'+self.B)
        qi_qj = self._qiqj
        if qi_qj == 'bu':
            Vij = ckm.get_ckm(self.pm)[0,2] # V_{ub} for b->u transitions
        if qi_qj == 'bc':
            Vij = ckm.get_ckm(self.pm)[1,2] # V_{cb} for b->c transitions
        if q2 <= ml**2:
            return 0
        
        p = GF*GF*(np.abs(Vij)**2)*(mB**5)/(192*(np.pi**3))
        return p

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M_SM(self, q2: float, mC: float = None) -> float:
        # From arxiv.org/pdf/1711.03110, eq. 31a
        # And https://scoap3-prod-backend.s3.cern.ch/media/files/79991/10.1088/1674-1137/ace821.pdf eq. 109
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.
        
        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC, mB)
        kA1 = ff["kA1"]
        kA2 = ff["kA2"]
        kA3 = ff["kA3"]
        kV  = ff["kV"]
        wsqm1 = (w*w - 1)
        wmr = w-r
        wmrsq = wmr*wmr

        gamma = (1./3.)*gamma0*(1/(mB*mC))*(r**3)*(wsqm1**1.5)*(q2hat - rhol)**2/(q2hat**3)*(
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

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2_SM(self, q2: float) -> float:
        return self.dGdq2M_SM(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def dGdq2M(self, q2: float, mC: float) -> float:
        gammaSM = self.dGdq2M_SM(q2, mC)
        gamma0 = self._rate_prefactor(q2)
        if gamma0 <= 0:
            return 0.

        wc = get_wceff_fccc(self.pm.get_wc(), self.pm, self._qiqj, self.lep, self.nu, self.scale, nf=5, withSM = False)

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

        mB = self.get_param("m_"+self.B)
        mC = self.get_param("m_"+self.M) if mC is None else mC
        r = mC/mB
        rhol = (self.get_param('m_'+self.lep)**2)/(mB**2)
        q2hat = q2/(mB**2)
        w = (mB**2 + mC**2 - q2) / (2 * mB * mC)
        if w < 1.:
            return 0.0
        
        ff = self.FF.get_ff_mmeson(q2, mC, mB)
        kA1 = ff["kA1"]
        kA2 = ff["kA2"]
        kA3 = ff["kA3"]
        kT1 = ff["kT1"]
        kT2 = ff["kT2"]
        kT3 = ff["kT3"]
        kV  = ff["kV"]
        kP  = ff["kP"]
        wsqm1 = (w**2 - 1)

        # This should follow https://arxiv.org/pdf/1711.03110 Eq. 32b, needs to be checked, high risk of typos
        nscale = (2./3.)*gamma0*(1/(mB*mC))*(r**3)*(wsqm1**1.5)*(q2hat - rhol)**2/(q2hat**3)
        gamma = gammaSM*(np.abs(1. + Cvl - Cvr)**2) + nscale*(
            6*(np.real(Cvr) + np.real(Cvr*np.conjugate(Cvl)))*kV*kV*wsqm1*q2hat*(2*q2hat + rhol)
            + 3*(np.abs(Csr - Csl)**2)*kP*kP*wsqm1*q2hat*q2hat
            + 6*kP*wsqm1*q2hat*np.sqrt(rhol)*(np.real(Csr-Csl) 
                    + np.real((Csr-Csl)*np.conjugate(Cvl)) 
                    - np.real((Csr-Csl)*np.conjugate(Cvr))
                )*(
                kA1 + kA2 (1 - r*w) + kA3*(w - r)
            )
            + 16*np.abs(Ct)**2 * (q2hat + 2*rhol)*(
                kT1*kT1*(w + 1)*(q2hat*(4*w + 1) + 6*r*wsqm1) 
                + kT2*kT2*(w - 1)*(q2hat*(4*w - 1) + 6*r*wsqm1)
                + kT3*wsqm1*q2hat*(
                    kT3*wsqm1 + 2*kT1*(w + 1) + 2*kT2*(w - 1)
                ) 
                - 4*kT1*kT2*wsqm1*(1 + r*w - 2*r*r)
            )
            + 12*np.sqrt(rhol)*q2hat*(
                wsqm1*(np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) - np.real(Ct*np.conjugate(Cvr)))*(
                    2*(kA2*r + kA3)*(kT1*(w + 1) + (w - 1)(kT2 + kT3*(1 + w)))
                )
                - 3*kV*wsqm1*(np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) + np.real(Ct*np.conjugate(Cvr)))*(
                    kT1*(1 + r) - kT2*(1 - r)
                )
                + kA1*(np.real(Ct) + np.real(Ct*np.conjugate(Cvl)) - np.real(Ct*np.conjugate(Cvr)))*(
                    kT1*(w + 1)*(3 + 2*w - 5*r) - kT2*(w - 1)*(3 - 2*w + 5*r) + 2*kT3*wsqm1*(w - r)
                )
            )
        )
        return gamma
    
    @fluctsettings(FluctType.NUMERIC)
    def dGdq2(self, q2: float) -> float:
        return self.dGdq2M(q2, None)

    @fluctsettings(FluctType.NUMERIC)
    def Gamma(self) -> float:
        return flavio.math.integrate.nintegrate(self.dGdq2, self.q2min, self.q2max)
        
