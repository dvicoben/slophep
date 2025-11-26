from slophep.Predictions.FormFactorsBaryonic import LbToLcFF
import check_utils as chk

import numpy as np
import eos


# Getting the EOS predictions:
slop_to_eos = {
    "T50" : "{}::f_long^T5(q2)",
    "T0"  : "{}::f_long^T(q2)",
    "A0"  : "{}::f_long^A(q2)",
    "V0"  : "{}::f_long^V(q2)",
    "Tp"  : "{}::f_perp^T(q2)",
    "T5p" : "{}::f_perp^T5(q2)",
    "Vp"  : "{}::f_perp^V(q2)",
    "Ap"  : "{}::f_perp^A(q2)",
    "At"  : "{}::f_time^A(q2)",
    "Vt"  : "{}::f_time^V(q2)"
}


def get_spectrum_eos(ffs: list, qsq, mode, ffscheme):
    options = eos.Options({
        'l' : "mu",
        'model' : 'WET',
        'form-factors' : ffscheme
    })
    ff = {}
    for iff in ffs:
        eosff = slop_to_eos[iff].format(mode)
        ffeval = [eos.Observable.make(eosff, eos.Parameters(), eos.Kinematics(q2=iq2), options).evaluate() for iq2 in qsq]
        ff[iff] = np.array(ffeval)
    return ff


# SLOP predictions
slopFF = LbToLcFF.DKMR()
q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")


# EOS predictions
eosFF_spectrum = get_spectrum_eos(["T50", "T0", "A0", "V0", "Tp", "T5p", "Vp", "Ap", "At", "Vt"], qsq, "Lambda_b->Lambda_c", "DKMR2017")

# Comparison plots
ff = ["T50", "T0", "A0", "V0", "Tp", "T5p", "Vp", "Ap", "At", "Vt"]
fflabel = [r"$f_\parallel^{T5}$", r"$f_\parallel^{T}$", r"$f_\parallel^{A}$", r"$f_\parallel^V$", 
    r"$f_\perp^T$", r"$f_\perp^{T5}$", r"$f_\perp^V$", r"$f_\perp^A$", r"$f_t^A$", r"$f_t^V$"]
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_eos_LbToLcFFDKMR_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_comparison_prediction(qsq, eosFF_spectrum[iff], "EOS")
    cplot.makeplot()
    cplot.savefig(savepath)