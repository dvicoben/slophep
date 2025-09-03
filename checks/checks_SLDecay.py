"""
Cross-check with SL_Decay, see
- https://cds.cern.ch/record/2313977/files/LHCb-INT-2018-015.pdf
- https://gitlab.cern.ch/scali/SL_Decay/-/tree/master?ref_type=heads

NOTE: for BGL
- SL_Decay has different defaults and internalparams, in particular, 
    - Different resonance masses for Blaschke factors,
      Additionally, for g, the heavies resonance is removed as considered to close to threshold
    - Different chim, chip, chimL
    - There is a bug in the outer functions, where pow(xx, 3/2) is used (similar issue with 5/2).
      This results in an int division of 3/2=1 and incorrect exponentation. 
      For the comparison this has been corrected. 
      A MR has been created in SL_Decay with the fix.
    - F2 has an additional factor of sqrt(r)/(1+r) where r=M(Charm meson)/M(B meson)
        - According to https://arxiv.org/pdf/1707.09509 F2 = (1+r)/(sqrt(r)) P1 so maybe
          this is meant to convert to P1?
"""
from slophep.Predictions.FormFactorsBToV import BdToDstFF

import numpy as np
import matplotlib.pyplot as plt

dat_BGL = np.loadtxt("checks/BdToDst_BGL_SLDecay.txt", float, skiprows=24)
dat_BGL = dat_BGL.T
q2 = dat_BGL[0]
BdToDstBGL_SLDecay = {
    "g"    : dat_BGL[2],
    "f"    : dat_BGL[3],
    "F1"   : dat_BGL[4],
    "F2"   : dat_BGL[5],
    "Blf"  : dat_BGL[6],
    "Blg"  : dat_BGL[7],
    "BlP1" : dat_BGL[8],
    "Phig" : dat_BGL[9],
    "Phif" : dat_BGL[10],
    "PhiF1" : dat_BGL[11],
    "PhiF2" : dat_BGL[12],
    "z" : dat_BGL[1]
}

SLDecay_BGL_ffpar = {
    "a0" : 0.0289,
    "a1" : 0.08,
    "a2" : -1,
    "b0" : 0.01224,
    "b1" : -0.052,
    "b2" : 1,
    "c1" : -0.007,
    "c2" : 0.089,
    "d0" : 0.0595,
    "d1" : -0.318
}

SLDecay_BGL_internalparams = {
    "Mb" : 5.27963,
    "chim"       : 3.894e-4,
    "chip"       : 5.131e-4,
    "chimL"      : 1.9421e-2,
    "BcStatesf"  : np.array([6.739, 6.750, 7.145, 7.150]),  # GeV
    # "BcStatesg"  : np.array([6.329, 6.920, 7.020, 7.280]),         # GeV
    "BcStatesg"  : np.array([6.329, 6.920, 7.020]),         # GeV
    "BcStatesP1" : np.array([6.275, 6.842, 7.250])          # GeV
}

BdToDstBGL = BdToDstFF.BGL()
BdToDstBGL.set_ff(**SLDecay_BGL_ffpar)
BdToDstBGL._internalparams.update(SLDecay_BGL_internalparams)
rC = BdToDstBGL.internalparams["Mc"]/BdToDstBGL.internalparams["Mb"]

def get_spectrum_BGL(qsq, ff):
    res = {}
    for iq2 in qsq:
        iff = ff.get_ff_gfF1F2_basis(iq2)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    res = {k : np.array(res[k]) for k in res}
    res["F2"] = res["F2"]*np.sqrt(rC)/(1+rC)
    return res


def make_comparison_plot(sloppred, otherpred, obsplot, qsq, prefix):
    for ipred in obsplot:
        fig, ax = plt.subplots(1, 1)
        ax.plot(qsq, sloppred[ipred], 'b-', label="SLOP")
        ax.plot(qsq, otherpred[ipred], 'r--', label="SL Decay")
        ax.set(xlabel = r"$q^2$", ylabel=ipred, title=prefix)
        ax.legend()
        plt.savefig(f"checks/checks_SLDecay_{prefix}_{ipred}.png", 
                    bbox_inches = 'tight',
                    dpi=100)


BdToDstBGL_ff = get_spectrum_BGL(q2, BdToDstBGL)

make_comparison_plot(BdToDstBGL_ff, BdToDstBGL_SLDecay, ["f", "g", "F1", "F2"], q2, "BdToDstFFBGL")
