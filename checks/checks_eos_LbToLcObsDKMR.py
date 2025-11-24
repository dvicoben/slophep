from slophep.Predictions.FormFactorsBaryonic import LbToLcFF
from slophep.Predictions.Observables import LbToLcEllNuPrediction
import check_utils as chk

import numpy as np
import eos


# Things for the EOS predictions:
slop_to_eos = {
    "afb_comb" : "{}::A_FB^c(q2)",
    "afb_had"  : "{}::A_FB^h(q2)",
    "afb_lep"  : "{}::A_FB^l(q2)",
    "f0"       : "{}::F_0(q2)",
    "dBRdq2"   : "{}::dBR/dq2",
    "1c"       : "{}::K_1c",
    "1cc"      : "{}::K_1cc",
    "1ss"      : "{}::K_1ss",
    "2c"       : "{}::K_2c",
    "2cc"      : "{}::K_2cc",
    "2ss"      : "{}::K_2ss",
    "3s"       : "{}::K_3s",
    "3sc"      : "{}::K_3sc",
    "4s"       : "{}::K_4s",
    "4sc"      : "{}::K_4sc",
}


def get_spectrum_eos(obs: list, qsq, mode, ffscheme):
    options = eos.Options({
        'l' : "mu",
        'model' : 'WET',
        'form-factors' : ffscheme
    })
    ff = {}
    for iff in obs:
        eosff = slop_to_eos[iff].format(mode)
        ffeval = [eos.Observable.make(eosff, eos.Parameters(), eos.Kinematics(q2=iq2), options).evaluate() for iq2 in qsq]
        ff[iff] = np.array(ffeval)
    return ff


def get_slop_spectrum(qsq, slopObs):
    slop_spectrum = {
        "afb_lep"  : np.array([slopObs.afb_lep(iq2) for iq2 in qsq]),
        "afb_had"  : np.array([slopObs.afb_had(iq2) for iq2 in qsq]),
        "afb_comb" : np.array([slopObs.afb_comb(iq2) for iq2 in qsq]),
        "f0"       : np.array([slopObs.f0(iq2) for iq2 in qsq]),
        "dBRdq2"   : np.array([slopObs.dBRdq2(iq2) for iq2 in qsq])
    }
    return slop_spectrum

# SLOP predictions
slopObs = LbToLcEllNuPrediction("mu", "mu", LbToLcFF.DKMR)
slopObs_aligned = LbToLcEllNuPrediction("mu", "mu", LbToLcFF.DKMR)
slopObs_aligned.set_alphaL(-0.78)
q2min, q2max = slopObs.q2min, slopObs.q2max
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slop_spectrum = get_slop_spectrum(qsq, slopObs)
slop_spectrum_aligned = get_slop_spectrum(qsq, slopObs_aligned)

# EOS predictions
eos_spectrum = get_spectrum_eos(["afb_lep", "afb_had", "afb_comb", "f0", "dBRdq2"], 
                                qsq, "Lambda_b->Lambda_clnu", "DKMR2017")

# Comparison plots
obs = ["afb_lep", "afb_had", "afb_comb", "f0", "dBRdq2"]
obslabel = [r"$A_{FB}^{\ell}$", r"$A_{FB}^{h}$", r"$A_{FB}^{c}$",
            r"$F_0$", r"$\mathrm{d}\mathcal{B}/\mathrm{d}q^2$"]
annotation = r"""Notes:
- SLOP (default) reproduces EOS 
  functionality using $\alpha_{-}^{\Lambda_c} = -0.786$
  from latest PDG value for $\Lambda_c^+ \to \Lambda\pi^+$
- For aligned we set $\alpha_{-}^{\Lambda_c} = -0.78$ 
  as in EOS
"""
for iobs, iobslabel in zip(obs, obslabel):
    savepath = f"checks/check_eos_LbToLcObsDKMR_{iobs}.png"
    cplot = chk.ComparisonPlot(iobslabel)
    cplot.add_slop_prediction(qsq, slop_spectrum[iobs], "SLOP (default)")
    cplot.add_slop_prediction(qsq, slop_spectrum_aligned[iobs], "SLOP (aligned)")
    cplot.add_comparison_prediction(qsq, eos_spectrum[iobs], "EOS")
    cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)