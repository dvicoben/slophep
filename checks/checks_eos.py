from slophep.Predictions.FormFactorsBToP import BToDFF
from slophep.Predictions.FormFactorsBToV import BToDstFF

import numpy as np
import matplotlib.pyplot as plt

# Getting the SLOP predictions for BSZ:
def get_spectrum(ff):
    q2max = (ff.internalparams["Mb"] - ff.internalparams["Mc"])**2
    qsq = np.linspace(0.012, q2max-1e-3, 100)
    
    res = {}
    for iq2 in qsq:
        iff = ff.get_ff(iq2)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return qsq, {k : np.array(res[k]) for k in res}


btod_bsz = BToDFF.BSZ()
btod_qsq, btod_bszff = get_spectrum(btod_bsz)

btodst_bsz = BToDstFF.BSZ()
btodst_qsq, btodst_bszff = get_spectrum(btodst_bsz)


# Getting the EOS predictions:
import eos
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


btod_bszeos = get_spectrum_eos(["f+", "f0", "fT"], btod_qsq, "B->D", "BSZ2015")
btodst_bszeos = get_spectrum_eos(["A0", "A1", "A12", "V", "T1", "T2", "T23"], btodst_qsq, "B->D^*", "BSZ2015")
# bstodsst_bszeos = get_spectrum_eos(["A0", "A1", "A12", "V", "T1", "T2", "T23"], btodst_qsq, "B_s->D_s^*", "BGJvD2019")


# Plotting them together
def make_comparison_plot(sloppred, eospred, qsq, prefix):
    for ipred in sloppred:
        fig, ax = plt.subplots(1, 1)
        ax.plot(qsq, sloppred[ipred], 'b-', label="SLOP BSZ")
        ax.plot(qsq, eospred[ipred], 'r--', label="EOS BSZ")
        ax.set(xlabel = r"$q^2$", ylabel=ipred)
        ax.legend()
        plt.savefig(f"checks/checks_eos_{prefix}{ipred}.png", 
                    bbox_inches = 'tight',
                    dpi=100)

make_comparison_plot(btod_bszff, btod_bszeos, btod_qsq, "BToDFF_")
make_comparison_plot(btodst_bszff, btodst_bszeos, btodst_qsq, "BToDstFF_")