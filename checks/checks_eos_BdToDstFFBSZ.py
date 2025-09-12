from slophep.Predictions.FormFactorsBToV import BdToDstFF
import check_utils as chk

import numpy as np
import eos


# Getting the EOS predictions:
slop_to_eos = {
    "f+"  : "{}::f_+(q2)",
    "f0"  : "{}::f_0(q2)",
    "fT"  : "{}::f_T(q2)",
    "A0"  : "{}::A_0(q2)",
    "A1"  : "{}::A_1(q2)",
    "A12" : "{}::A_12(q2)",
    "V"   : "{}::V(q2)",
    "T1"  : "{}::T_1(q2)",
    "T2"  : "{}::T_2(q2)",
    "T23" : "{}::T_23(q2)",
}

def get_spectrum_eos(ffs: list, qsq, mode, ffscheme):
    options = eos.Options({
        'l' : "mu",
        'q' : 'd',
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
slopFF = BdToDstFF.BSZ()
q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

# EOS predictions
eosFF_spectrum = get_spectrum_eos(["A0", "A1", "A12", "V", "T1", "T2", "T23"], qsq, "B->D^*", "BSZ2015")

# Comparison plots
ff = ["A0", "A1", "A12", "V", "T1", "T2", "T23"]
fflabel = [r"$A_0$", r"$A_1$", r"$A_{12}$", r"$V$", r"$T_1$", r"$T_2$", r"$T_{23}$"]
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_eos_BdToDstFFBSZ_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_comparison_prediction(qsq, eosFF_spectrum[iff], "EOS")
    cplot.makeplot()
    cplot.savefig(savepath)