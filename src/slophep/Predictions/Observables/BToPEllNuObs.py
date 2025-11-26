from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

from slophep.Predictions.FormFactorsBToP import FormFactorBToP
from slophep.Predictions.Observables import ObservableBase
import slophep.Predictions.Math.BToPMathTools as mt

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
        """Alias for get\_angularcoeff

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficents a, b, c
        """
        return self.get_angularcoeff(q2)

    def J(self, q2: float) -> float:
        """Calculate rate normalised angular coefficients, flavio is used for these calculations, method essentially follows
        flavio.physics.bdecays.bplnu._get_angularcoeff from https://flav-io.github.io/apidoc/flavio/physics/bdecays/bplnu.m.html

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            Dictionary of coefficients a, b, c
        """
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
    
    def dGdq2_bin(self, q2min: float, q2max: float) -> float:
        """Caclulate binned q2 distriution

        Parameters
        ----------
        q2min : float
        q2max : float

        Returns
        -------
        float
            dGamma/dq2 (up to normalisation) integrated over the bin
        """
        J = self.dJ_bin(q2min, q2max)
        return 2 * (J['a'] + J['c']/3.)

    def dGdq2_hist(self, q2_bins: int | list):
        """Create 1D histogram of dG/dq2

        Parameters
        ----------
        q2_bins : int | list
            Binning in q2, specify either number of bins or bin edges. If an int is not provided it is
            assumed an iterable for bin edges has been provided.
        
        Returns
        -------
        h : list[float]
            PDF histogram
        q2_edges : list[float]
            List of bin edges used
        """
        q2_edges = q2_bins if type(q2_bins) != int else np.linspace(self.q2min, self.q2max, q2_bins+1, endpoint=True)
        h = np.zeros(len(q2_edges)-1)
        for iq2 in range(len(q2_edges)-1):
            idG = self.dGdq2_bin(q2_edges[iq2], q2_edges[iq2+1])
            h[iq2] = idG
        return h, q2_edges

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
            Dictionary of observables, integrated over q2min, q2max
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
            Dictionary of observables, integrated over q2min, q2max
        """
        dJ = self.dJ_bin(q2min, q2max)
        norm = 2 * (dJ['a'] + dJ['c']/3.)
        # den = flavio.math.integrate.nintegrate(self.dGdq2, q2min, q2max)
        # return {iobs : num[iobs]/den for iobs in self.obslist}
        return {iobs : dJ[iobs]/norm for iobs in self.obslist}

    def dJ_q2int(self) -> dict:
        """Calculate q2-integrated observable

        Returns
        -------
        dict
            Dictionary of observables
        """
        return self.dJ_bin(self.q2min, self.q2max)

    def J_q2int(self) -> dict:
        """Calculate rate-normalised q2-integrated observable

        Returns
        -------
        dict
            Dictionary of observables
        """
        dJ = self.dJ_q2int()
        norm = 2 * (dJ['a'] + dJ['c']/3.)
        return {iobs : dJ[iobs]/norm for iobs in self.obslist}

    def afb(self, q2: float) -> float:
        """Calculate afb"""
        return self.J(q2)["b"]
    
    def afb_bin(self, q2min: float, q2max: float) -> float:
        """Calculate binned AFB"""
        j = self.J_bin(q2min, q2max)
        return j["b"]
    
    def PDF(self, q2: float, ctl: float) -> float:
        """Evaluate PDF (up to normalisation) at phase-space point

        Parameters
        ----------
        q2 : float
        ctl : float

        Returns
        -------
        float
            pdf (up to normalisation)
        """
        j = self.dJ(q2)
        return mt.angularPDF(ctl, j)
    
    def PDF_norm(self, q2: float, ctl: float) -> float:
        """Evaluate PDF (up to normalisation) at phase-space point, using rate-normalised observables

        Parameters
        ----------
        q2 : float
        ctl : float

        Returns
        -------
        float
            pdf (up to normalisation)
        """
        j = self.J(q2)
        return mt.angularPDF(ctl, j)

    def PDF_bin(self, q2_min: float, q2_max: float,
                 ctl_min: float, ctl_max: float) -> float:
        """Evaluate pdf integrated over phase-space bin

        Parameters
        ----------
        q2_min : float
        q2_max : float
        ctl_min : float
        ctl_max : float

        Returns
        -------
        float
            PDF in phase-space bin
        """
        j = self.dJ_bin(q2_min, q2_max)
        return mt.angularPDF_binned(ctl_min, ctl_max, j)
    
    def PDF_angular_int(self, ctl_min: float, ctl_max: float) -> dict[float]:
        """Evaluate angular terms of PDF integrated over angular bin
        
        Parameters
        ----------
        ctl_min : float
        ctl_max : float

        Returns
        -------
        dict
            Integrated angular term corresponding to each observable
        """
        return mt.angular_integrals(ctl_min, ctl_max)
    
    def PDF_hist_angular_int(self, ctl_bins: int | list):
        """Compute angular integrals for provided binning scheme

        Parameters
        ----------
        ctl_bins : int | list
            Binning in ctl, specify either number of bins or bin edges

        Returns
        -------
        dict
            Angular integrals for angular bin, output is in the form dict[ctx bin][ctl bin][chi bin] = [angular integrals]
            where [angular integrals] is a list of the integrated angular term corresponding to each observable in the order
            stored in dict["order"]. This is also returned as an array in dict["asarray"] to use for easier multiplication with
            numpy.
        """
        ctl_edges = ctl_bins if type(ctl_bins) != int else np.linspace(-1.0, 1.0, ctl_bins+1, endpoint=True)
        h_angintegrals = {jctl : {} for jctl in range(len(ctl_edges)-1)}
        angint_array = np.zeros(len(ctl_edges)-1, 3)
        for ictl in range(len(ctl_edges)-1):
            h_ang_d = self.PDF_angular_int(ctl_edges[ictl], ctl_edges[ictl+1])
            h_angintegrals[ictl] = np.array([h_ang_d[k] for k in self.obslist])
            angint_array[ictl] = np.array([h_ang_d[k] for k in self.obslist])
        h_angintegrals["order"] = self.obslist.copy()
        h_angintegrals["asarray"] = angint_array
        return h_angintegrals

    def PDF_hist(self, q2_bins: int | list, ctl_bins: int | list):
        """Create 4D histogram of PDF. This computes observables integrated over the bin and angular integrals.

        Parameters
        ----------
        q2_bins : int | list
            Binning in q2, specify either number of bins or bin edges
        ctl_bins : int | list
            Binning in ctl, specify either number of bins or bin edges

        Returns
        -------
        h : np.ndarray
            PDF histogram
        bin_edges : list
            List of bin edges used, in form [q2_edges, ctd_edges, ctl_edges, chi_edges]
        h_angint : dict
            Angular integrals for angular bin
        j_bins : dict
            Observables computed to generate h
        """
        q2max = self.q2max
        q2min = self.q2min
        
        q2_edges = q2_bins if type(q2_bins) != int else np.linspace(q2min, q2max, q2_bins+1, endpoint=True)
        ctl_edges = ctl_bins if type(ctl_bins) != int else np.linspace(-1.0, 1.0, ctl_bins+1, endpoint=True)

        h = np.zeros((len(q2_edges)-1, len(ctl_edges)-1))
        # Cache angular integrals to avoid recomputation for q2 angular bin
        h_angintegrals = self.PDF_hist_angular_int(ctl_edges)
        j_bins = {}
        for iq2 in range(len(q2_edges)-1):
            for ictl in range(len(ctl_edges)-1):
                dJ = self.dJ_bin(q2_edges[iq2], q2_edges[iq2+1])
                j_bins[iq2] = dJ
                J_vec = np.array([dJ[iobs] for iobs in h_angintegrals["order"]])
                h[iq2][ictl] = 9/(32*np.pi)*np.dot(J_vec, h_angintegrals[ictl])
        return h, [q2_edges, ctl_edges], h_angintegrals, j_bins
    
    def plot_obs_prediction(self, obs: str, 
                            q2min: float = None, q2max: float = None, 
                            label: str = None,
                            plot: tuple[plt.Figure, plt.Axes] = None) -> tuple[plt.Figure, plt.Axes]:
        """Plot prediction for a particular observable

        Parameters
        ----------
        obs : str
            Desired observable, available are ["a", "b", "c", "AFB"]
        q2min : float, optional
            Min. value of q2 to plot, by default None which sets to physical minimum
        q2max : float, optional
            Max. value of q2 to plot, by default None which sets to physical maximum
        label : str, optional
            Label for legend, by default None
        plot : tuple[plt.Figure, plt.Axes], optional
            Figure and axes to plot on, by default None. This allows to plot multiple observables in same
            axis or to change FFs/WCs and re-plot on same axes

        Returns
        -------
        fig : plt.Figure
            The matplotlib figure, show with fig.show(), update drawing with fig.canvas.draw()
        ax : plt.Axes
            Matplotlib axes, add legend with ax.legend()

        Raises
        ------
        ValueError
            For unavailable observable
        """
        q2min = q2min if q2min else self.q2min+1e-6
        q2max = q2max if q2max else self.q2max-1e-6
        q2 = np.linspace(q2min, q2max, 150, endpoint=True)
        obs_vals = []
        # Get the observable in q2
        if obs == "AFB":
            obs_vals = np.array([self.afb(iq2) for iq2 in q2])
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