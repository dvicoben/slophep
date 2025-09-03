"""
Comparison in BLPR with flavio

This requires a few workaronds. The calculation of BLPR in SLOP is lifted from HAMMER.
- flavio uses different quark masses, lambdabar and has subsubleading order contributions.
- flavio also calculates the leading IW function Xi completely differently.

For comparison:
- We look at the hatted h functions as aligning the leading IW Xi is not easy
- HAMMER has an ebReb and ecRec correction that is not present in flavio. We set it to zero in SLOP for this comparison.
- All parameters for flavio are set to the ones in SLOP
- All subsubleading contributions in flavio are set to 0.

NOTE: This still does not completely align things, as L_i are calculated differently in flavio 
(with w-1 determined from an expansion z around z=0) so exact match is not expected, and there
may be small disagreements in the very low q2
"""

from slophep.Predictions.FormFactorsBToV import BdToDstFF
from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt
import flavio
import check_utils as chk

setPlotParams()

# Getting the SLOP predictions for BLPR:
slopFF = BdToDstFF.BLPR()
slopFF.internalparams.update({
    "ebReb" : 1.0,
    "ecRec" : 1.0
})
q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_hhat")

def align_pars_for_flavio(par: dict, ff):
    ffpar = ff.ffpar
    par = par.copy()
    par["chi_2(1)"] = ffpar["Chi21"]
    par["chi_2p(1)"] = ffpar["Chi2p"]
    par["chi_2pp(1)"] = 0.0
    par["chi_3p(1)"] = ffpar["Chi3p"]
    par["chi_3pp(1)"] = 0.0
    par["eta(1)"] = ffpar["Eta1"]
    par["etap(1)"] = ffpar["Etap"]
    par["etapp(1)"] = 0.0
    for i in range(1, 7):
        par[f"CLN l_{i}(1)"] = 0.0
        par[f"CLN lp_{i}(1)"] = 0.0
    par["__ash"] = ff.internalparams["ash"]
    la = ff.internalparams["la"]
    mb = ff.internalparams["mb"]
    eb = la/(2*mb)
    mc = mb - ff.internalparams["delta_mbc"]
    ec = la/(2*mc)
    par["__mc"] = mc
    par["__epsc"] = ec
    par["__mb"] = mb
    par["__epsb"] = eb
    return par

par = flavio.default_parameters.get_central_all()
par = align_pars_for_flavio(par, slopFF)

from flavio.physics.bdecays.formfactors import hqet
from flavio.physics.bdecays.formfactors import common

process_dict = {
    'B->D'  : {'B': 'B+', 'q': 'b->c', 'P': 'D0'},
    'B->D2' : {'B': 'B0', 'q': 'b->c', 'P': 'D+'},
    'B->D*' : {'B': 'B0', 'V': 'D*+', 'q': 'b->c'}
}
# We lift this here from flavio and adapt it for our purposes
def blpr_btov_flavio(process, q2, par, scale, order_z=3, order_z_slp=2, order_z_sslp=1):
    r"""Central value of $B\to V$ form factors in the lattice convention
    CLN parametrization.

    See arXiv:1703.05330.
    """
    pd = process_dict[process]
    mB = par['m_' + pd['B']]
    mV = par['m_' + pd['V']]
    w = max((mB**2 + mV**2 - q2) / (2 * mB * mV), 1)
    # We also need to match the quark masses which differ in flavio
    # phqet = hqet.get_hqet_parameters(par)
    ash = par["__ash"]
    epsc = par["__epsc"]
    epsb = par["__epsb"]
    zc = par["__mc"]/par["__mb"]
    # eq. (22) of arXiv:0809.0222
    CV1 = hqet.CV1(w, zc)
    CA1 = hqet.CA1(w, zc)
    CA2 = hqet.CA2(w, zc)
    CA3 = hqet.CA3(w, zc)
    CT1 = hqet.CT1(w, zc)
    CT2 = hqet.CT2(w, zc)
    CT3 = hqet.CT3(w, zc)
    z = common.z(mB, mV, q2, t0='tm')
    rho2 = par['CLN rho2_xi']
    c = par['CLN c_xi']
    xi3 = par['CLN xi3']
    xi = hqet.xi(z, rho2, c, xi3, order_z=order_z)
    L = hqet.Lz(par, w, z, order_z=order_z_slp)
    ell = hqet.ell(par, z, order_z=order_z_sslp)
    h = {}
    h['V'] = (1 + ash * CV1
              + epsc * (L[2] - L[5])
              + epsb * (L[1] - L[4])
              + epsc**2 * (ell[2] - ell[5]))
    h['A1'] = (1 + ash * CA1
               + epsc * (L[2] - L[5] * (w - 1)/(w + 1))
               + epsb * (L[1] - L[4] * (w - 1)/(w + 1))
               + epsc**2 * (ell[2] - ell[5] * (w - 1)/(w + 1)))
    h['A2'] = (ash * CA2
               + epsc * (L[3] + L[6])
               + epsc**2 * (ell[3] + ell[6]))
    h['A3'] = (1 + ash * (CA1 + CA3)
               + epsc * (L[2] - L[3] + L[6] - L[5])
               + epsb * (L[1] - L[4])
               + epsc**2 * (ell[2] - ell[3] + ell[6] - ell[5]))
    h['T1'] = (1 + ash * (CT1 + (w - 1)/2 * (CT2 - CT3))
               + epsc * L[2]
               + epsb * L[1]
               + epsc**2 * ell[2])
    h['T2'] = (ash * (w + 1)/2 * (CT2 + CT3)
               + epsc * L[5]
               - epsb * L[4]
               + epsc**2 * ell[5])
    h['T3'] = (ash * CT2
               + epsc * (L[6] - L[3])
               + epsc**2 * (ell[6] - ell[3]))
    return h

def get_flavio_spectrum(fffunc, process, par, qsq):
    res = {}
    for iq2 in qsq:
        iff = fffunc(process, iq2, par, 4.8)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return {k : np.array(res[k]) for k in res}

flavioFF_spectrum = get_flavio_spectrum(blpr_btov_flavio, "B->D*", par, qsq)


chk.make_comparison_plot(slopFF_spectrum, flavioFF_spectrum, qsq, "flavio",
    ["A1", "A2", "A3", "T1", "T2", "T3", "V"],
    [r"$\hat{h}_{A1}$", r"$\hat{h}_{A2}$", r"$\hat{h}_{A3}$", 
     r"$\hat{h}_{T1}$", r"$\hat{h}_{T2}$", r"$\hat{h}_{T3}$", r"$\hat{h}_{V}$"],
    "BdToDstFFBLPR", "checks/check_flavio_{}_{}.png")
