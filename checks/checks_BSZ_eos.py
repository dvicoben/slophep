from slophep.Predictions.FormFactorsBToP import BuToDFF, BdToPiFF
from slophep.Predictions.FormFactorsBToV import BdToDstFF

import numpy as np
import matplotlib.pyplot as plt
import eos

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




# Plotting them together
def make_comparison_plot(sloppred, eospred, qsq, prefix):
    for ipred in sloppred:
        fig, ax = plt.subplots(1, 1)
        ax.plot(qsq, sloppred[ipred], 'b-', label="SLOP BSZ")
        ax.plot(qsq, eospred[ipred], 'r--', label="EOS BSZ")
        ax.set(xlabel = r"$q^2$", ylabel=ipred, title=prefix)
        ax.legend()
        plt.savefig(f"checks/checks_eos_{prefix}_{ipred}.png", 
                    bbox_inches = 'tight',
                    dpi=100)


if __name__ == "__main__":
    # SLOP predictions
    btodst_bsz = BdToDstFF.BSZ()
    btodst_qsq, btodst_bszff = get_spectrum(btodst_bsz)
    btod_bsz = BuToDFF.BSZ()
    btod_qsq, btod_bszff = get_spectrum(btod_bsz)
    btopi_bsz = BdToPiFF.BSZ()
    btopi_qsq, btopi_bszff = get_spectrum(btopi_bsz)
    
    # EOS predictions
    btodst_bszeos = get_spectrum_eos(["A0", "A1", "A12", "V", "T1", "T2", "T23"], btodst_qsq, "B->D^*", "BSZ2015")
    btod_bszeos = get_spectrum_eos(["f+", "f0", "fT"], btod_qsq, "B->D", "BSZ2015")
    btopi_bszeos = get_spectrum_eos(["f+", "f0", "fT"], btopi_qsq, "B->pi", "BSZ2015")

    # Comparison plots
    make_comparison_plot(btodst_bszff, btodst_bszeos, btodst_qsq, "BdToDstFFBSZ")
    make_comparison_plot(btod_bszff, btod_bszeos, btod_qsq, "BuToD0FFBSZ")
    make_comparison_plot(btopi_bszff, btopi_bszeos, btopi_qsq, "BdToPiFFBSZ")