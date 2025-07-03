from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

from bd2dstlnu.Predictions.BToDstFFBase import FormFactor
from bd2dstlnu.Predictions import BToDstMathTools as mt

import flavio
from flavio.physics.running import running
from flavio.physics.bdecays.wilsoncoefficients import get_wceff_fccc
import flavio.physics.bdecays.bvlnu as bvlnu
from flavio.physics.bdecays import angular



class BToDstEllNuPrediction:
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactor,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        self._B: str = "B0"
        self._V: str = "D*+"
        self._qiqj: str = "bc"
        self._lep: str = lep
        self._nu: str = nu
        self._wc_obj: flavio.WilsonCoefficients = flavio.WilsonCoefficients()

        self._par: dict = flavio.default_parameters.get_central_all()
        if type(par) == dict:
            self._par = par
        self._scale: float = scale

        self._FF: FormFactor = FF(par, scale, *ffargs)
        self.obslist = ['1s', '1c', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]

    @property
    def scale(self) -> float: 
        """Renorm scale"""
        return self._scale
    @property
    def par(self) -> dict: 
        """Dictionary of parameters, defaults to `flavio.default_parameters.get_central_all()`"""
        return self._par
    @property
    def lep(self) -> str: 
        """Lepton flavour (mu/e/tau)"""
        return self._lep
    @property
    def nu(self) -> str: 
        """Neutrino flavour (nu/e/tau)"""
        return self._nu
    @property
    def wc_obj(self) -> flavio.WilsonCoefficients: 
        """Wilson coefficient object"""
        return self._wc_obj
    @property
    def FF(self) -> FormFactor: 
        """Form factors"""
        return self._FF

    def set_wc(self, wc_dict: dict, eft: str = 'WET', basis: str = 'flavio'):
        """Set the wilson coefficients

        Parameters
        ----------
        wc_dict : dict
            Dictionary of wilson coefficients
        eft : str, optional
            EFT, by default 'WET'
        basis : str, optional
            WC basis (see https://wcxf.github.io/bases.html), by default 'flavio'
        """
        self._wc_obj.set_initial(wc_dict, self.scale, eft, basis)

    def set_ff(self, ffparams: dict):
        """Set form factor parameters

        Parameters
        ----------
        ffparams : dict
            Dictionary of form factor parameters, names should match FF.ffpars dictionary in the particular 
            scheme being used
        """
        self._FF.set_ff(**ffparams)

    def get_angularcoeff(self, q2: float) -> dict:
        """Calculate angular coefficients, flavio is used for these calculations, method essentially follows
        flavio.physics.bdecays.bvlnu._get_angularcoeff from https://flav-io.github.io/apidoc/flavio/physics/bdecays/bvlnu.m.html

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents J_i
        """
        mb = running.get_mb(self.par, self.scale)
        wc = get_wceff_fccc(self.wc_obj, self.par, self._qiqj, self.lep, self.nu, mb, self.scale, nf=5)
        if self.lep != self.nu and all(C == 0 for C in wc.values()):
            # if all WCs vanish, so does the AC!
            return {k: 0 for k in
                    ['1s', '1c', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]}
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self._B]
        mV = self.par['m_'+self._V]
        mlight = running.get_mc(self.par, self.scale) # this is needed for scalar contributions
        N = bvlnu.prefactor(q2, self.par, self._B, self._V, self.lep)
        ff = self.FF.get_ff(q2)
        h = angular.helicity_amps_v(q2, mB, mV, mb, mlight, ml, 0, ff, wc, N)
        J = angular.angularcoeffs_general_v(h, q2, mB, mV, mb, mlight, ml, 0)
        return J

    def get_norm_coeff(self, q2: float) -> dict:
        """Calculate rate normalised angular coefficients, flavio is used for these calculations, method essentially follows
        flavio.physics.bdecays.bvlnu._get_angularcoeff from https://flav-io.github.io/apidoc/flavio/physics/bdecays/bvlnu.m.html

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficients J_i
        """
        J = self.get_angularcoeff(q2)
        norm =  3/4. * (2 * J['1s'] + J['1c']) - 1/4. * (2 * J['2s'] + J['2c'])
        return {k : J[k]/norm for k in J}
    
    def dJ(self, q2: float) -> dict:
        """Alias for get\_angularcoeff

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents J_i
        """
        return self.get_angularcoeff(q2)

    def J(self, q2: float) -> dict:
        """Alias for get\_norm\_coeff

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficients J_i
        """
        return self.get_norm_coeff(q2)
    
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
        return 3/4. * (2 * J['1s'] + J['1c']) - 1/4. * (2 * J['2s'] + J['2c'])

    def _obsq2Bin(self, obs: str | int, q2min: float, q2max: float) -> float:
        def evalObs(q2):
            return self.dJ(q2)[obs]
        return flavio.math.integrate.nintegrate(evalObs, q2min, q2max)

    def dJ_bin(self, q2min: float, q2max: float) -> dict:
        """Calculate binned angular observable

        Parameters
        ----------
        q2min : float
        q2max : float

        Returns
        -------
        dict
            Dictionary of observable <J_i>, integrated over q2min, q2max
        """
        return {iobs: self._obsq2Bin(iobs, q2min, q2max) for iobs in self.obslist}

    def J_bin(self, q2min: float, q2max: float) -> dict:
        """Calculate rate normalised binned angular observable

        Parameters
        ----------
        q2min : float
        q2max : float

        Returns
        -------
        dict
            Dictionary of observable <J_i>, integrated over q2min, q2max
        """
        J = self.dJ_bin(q2min, q2max)
        norm = 3/4. * (2 * J['1s'] + J['1c']) - 1/4. * (2 * J['2s'] + J['2c'])
        # den = flavio.math.integrate.nintegrate(self.dGdq2, q2min, q2max)
        # return {iobs : num[iobs]/den for iobs in self.obslist}
        return {iobs : J[iobs]/norm for iobs in self.obslist}

    def dJ_q2int(self) -> dict:
        """Calculate q2-integrated observable

        Returns
        -------
        dict
            Dictionary of observable <J_i>
        """
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self._B]
        mV = self.par['m_'+self._V]
        q2max = (mB-mV)**2
        q2min = ml**2
        return self.dJ_bin(q2min, q2max)
    
    def J_q2int(self) -> dict:
        """Calculate rate-normalised q2-integrated observable

        Returns
        -------
        dict
            Dictionary of observable <J_i>
        """
        J = self.dJ_q2int()
        norm = 3/4. * (2 * J['1s'] + J['1c']) - 1/4. * (2 * J['2s'] + J['2c'])
        return {iobs : J[iobs]/norm for iobs in self.obslist}
    
    def PDF(self, q2: float, ctx: float, ctl: float, chi: float) -> float:
        """Evaluate 4D PDF (up to normalisation) at phase-space point

        Parameters
        ----------
        q2 : float
        ctx : float
        ctl : float
        chi : float

        Returns
        -------
        float
            pdf (up to normalisation)
        """
        j = self.dJ(q2)
        return mt.angularPDF(ctx, ctl, chi, j)
    
    def PDF_norm(self, q2: float, ctx: float, ctl: float, chi: float) -> float:
        """Evaluate 4D PDF at phase-space point, using rate-normalised observables

        Parameters
        ----------
        q2 : float
        ctx : float
        ctl : float
        chi : float

        Returns
        -------
        float
            pdf (up to normalisation)"""
        j = self.J(q2)
        return mt.angularPDF(ctx, ctl, chi, j)
    
    def PDF_bin(self, q2_min: float, q2_max: float,
                ctx_min: float, ctx_max: float,
                ctl_min: float, ctl_max: float,
                chi_min: float, chi_max: float) -> float:
        """Evaluate 4D PDF integrated over phase-space bin

        Parameters
        ----------
        q2_min : float
        q2_max : float
        ctx_min : float
        ctx_max : float
        ctl_min : float
        ctl_max : float
        chi_min : float
        chi_max : float

        Returns
        -------
        float
            PDF in phase-space bin
        """
        j = self.dJ_bin(q2_min, q2_max)
        
        return mt.angularPDF_binned(ctx_min, ctx_max, ctl_min, ctl_max, chi_min, chi_max, j)
    
    def PDF_norm_bin(self, q2_min: float, q2_max: float,
                     ctx_min: float, ctx_max: float,
                     ctl_min: float, ctl_max: float,
                     chi_min: float, chi_max: float) -> float:
        """Evaluate 4D PDF integrated over phase-space bin, using rate-normalised angular observables

        Parameters
        ----------
        q2_min : float
        q2_max : float
        ctx_min : float
        ctx_max : float
        ctl_min : float
        ctl_max : float
        chi_min : float
        chi_max : float

        Returns
        -------
        float
            PDF in phase-space bin
        """
        j = self.J_bin(q2_min, q2_max)
        return mt.angularPDF_binned(ctx_min, ctx_max, ctl_min, ctl_max, chi_min, chi_max, j)

    def PDF_angular_int(self, ctx_min: float, ctx_max: float,
                     ctl_min: float, ctl_max: float,
                     chi_min: float, chi_max: float) -> dict[float]:
        """Evaluate angular terms of PDF integrated over angular bin
        
        Parameters
        ----------
        ctx_min : float
        ctx_max : float
        ctl_min : float
        ctl_max : float
        chi_min : float
        chi_max : float

        Returns
        -------
        dict
            Integrated angular term corresponding to each observable
        """
        return mt.angular_integrals(ctx_min, ctx_max, ctl_min, ctl_max, chi_min, chi_max)

    def PDF_hist(self, q2_bins: int | list, ctx_bins: int | list, ctl_bins: int | list, chi_bins: int | list):
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self._B]
        mV = self.par['m_'+self._V]
        q2max = (mB-mV)**2
        q2min = ml**2
        
        q2_edges = q2_bins if type(q2_bins) != int else np.linspace(q2min, q2max, q2_bins+1, endpoint=True)
        ctx_edges = ctx_bins if type(ctx_bins) != int else np.linspace(-1.0, 1.0, ctx_bins+1, endpoint=True)
        ctl_edges = ctl_bins if type(ctl_bins) != int else np.linspace(-1.0, 1.0, ctl_bins+1, endpoint=True)
        chi_edges = chi_bins if type(chi_bins) != int else np.linspace(-np.pi, np.pi, chi_bins+1, endpoint=True)

        h = np.zeros((len(q2_edges)-1, len(ctx_edges)-1, len(ctl_edges)-1, len(chi_edges)-1))
        # Cache angular integrals to avoid recomputation for q2 angular bin
        h_angintegrals = self.PDF_hist_angular_int(ctx_edges, ctl_edges, chi_edges)
        for iq2 in range(len(q2_edges)-1):
            for ictx in range(len(ctx_edges)-1):
                for ictl in range(len(ctl_edges)-1):
                    for ichi in range(len(chi_edges)-1):
                        dJ = self.dJ_bin(q2_edges[iq2], q2_edges[iq2+1])
                        J_vec = np.array([dJ[iobs] for iobs in h_angintegrals["order"]])
                        h[iq2][ictx][ictl][ichi] = 9/(32*np.pi)*np.dot(J_vec, h_angintegrals[ictx][ictl][ichi])
                        # This is inefficient as recomputes angular integrals for every calculation
                        # h[iq2][ictx][ictl][ichi] = self.PDF_bin(q2_edges[iq2], q2_edges[iq2+1],
                        #                                         ctx_edges[ictx], ctx_edges[ictx+1],
                        #                                         ctl_edges[ictl], ctl_edges[ictl+1],
                        #                                         chi_edges[ichi], chi_edges[ichi+1])
                        
        return h, [q2_edges, ctx_edges, ctl_edges, chi_edges], h_angintegrals

    def PDF_hist_angular_int(self, ctx_bins: int | list, ctl_bins: int | list, chi_bins: int | list):
        ctx_edges = ctx_bins if type(ctx_bins) != int else np.linspace(-1.0, 1.0, ctx_bins+1, endpoint=True)
        ctl_edges = ctl_bins if type(ctl_bins) != int else np.linspace(-1.0, 1.0, ctl_bins+1, endpoint=True)
        chi_edges = chi_bins if type(chi_bins) != int else np.linspace(-np.pi, np.pi, chi_bins+1, endpoint=True)
        h_angintegrals = {jctx : {jctl : {} for jctl in range(len(ctl_edges)-1)} for jctx in range(len(ctx_edges)-1)}
        for ictx in range(len(ctx_edges)-1):
            for ictl in range(len(ctl_edges)-1):
                for ichi in range(len(chi_edges)-1):
                    h_ang_d = self.PDF_angular_int(ctx_edges[ictx], ctx_edges[ictx+1],
                                                        ctl_edges[ictl], ctl_edges[ictl+1],
                                                        chi_edges[ichi], chi_edges[ichi+1])
                    h_angintegrals[ictx][ictl][ichi] = np.array([h_ang_d[k] for k in self.obslist])
        h_angintegrals["order"] = self.obslist.copy()
        return h_angintegrals


    def afb(self, q2: float) -> float:
        """ Calculate Afb
            NOTE: angular convention for 6s and 6c means this may differ by a sign """
        j = self.dJ(q2)
        num = (3.0/8)*(j["6c"] + 2*j["6s"])
        denom = 3/4. * (2 * j['1s'] + j['1c']) - 1/4. * (2 * j['2s'] + j['2c'])
        return num/denom
    
    def afb_bin(self, q2min: float, q2max: float) -> float:
        """Calculate binned AFB"""
        j = self.dJ_bin(q2min, q2max)
        num = (3.0/8)*(j["6c"] + 2*j["6s"])
        denom = 3/4. * (2 * j['1s'] + j['1c']) - 1/4. * (2 * j['2s'] + j['2c'])
        return num/denom

    def fl(self, q2: float) -> float:
        """ Calculate FL """
        j = self.dJ(q2)
        num = 3*j["1c"] - j["2c"]
        denom = 3*(j["1c"] + 2*j["1s"]) - (j["2c"] + 2*j["2s"])
        return num/denom
    
    def fl_bin(self, q2min: float, q2max: float) -> float:
        """Calculate binned FL"""
        j = self.dJ_bin(q2min, q2max)
        num = 3*j["1c"] - j["2c"]
        denom = 3*(j["1c"] + 2*j["1s"]) - (j["2c"] + 2*j["2s"])
        return num/denom
    
    def uniang_obs(self, q2: float) -> dict[float]:
        """Calculate uniangular observables (FL, AFB, Flt, J3, J9)"""
        j = self.J(q2)
        return mt.calc_unaing_obs(j)
    
    def binuniang_obs(self, q2min: float, q2max: float) -> dict[float]:
        """Calculate binned uniangular observables (FL, AFB, Flt, J3, J9)"""
        j = self.J_bin(q2min, q2max)
        return mt.calc_unaing_obs(j)

    def plot_obs_prediction(self, obs: str, 
                            q2min: float = None, q2max: float = None, 
                            label: str = None,
                            plot = None) -> tuple[plt.Figure, plt.Axes]:
        ml = self.par['m_'+self.lep]
        mB = self.par['m_'+self._B]
        mV = self.par['m_'+self._V]
        
        q2max = q2max if q2max else (mB-mV)**2
        q2min = q2min if q2min else ml**2
        q2 = np.linspace(q2min, q2max, 150, endpoint=True)
        obs_vals = []
        # Get the observable in q2
        if obs == "FL":
            obs_vals = np.array([self.fl(iq2) for iq2 in q2])
        elif obs == "AFB":
            obs_vals = np.array([self.afb(iq2) for iq2 in q2])
        elif obs == "FLt":
            obs_vals = np.array([self.uniang_obs(iq2)["FLt"] for iq2 in q2])
        else:
            if obs not in self.obslist:
                raise ValueError(f"{obs} is not a valid observable")
            obs_vals = np.array([self.J(iq2)[obs] for iq2 in q2])

        fig, ax = None, None
        if type(plot) == type(None):
            fig, ax = plt.subplots(1, 1)
            ax.set(xlabel=r"$q^2$ [GeV$^2$]", xlim=(q2min, q2max))
        else:
            fig, ax = plot

        ax.plot(q2, obs_vals, label=label if label else f"{obs} {self.FF.name}")
        return fig, ax

