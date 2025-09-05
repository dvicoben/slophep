"""
Cross-check with SL_Decay, see
- https://cds.cern.ch/record/2313977/files/LHCb-INT-2018-015.pdf
- https://gitlab.cern.ch/scali/SL_Decay/-/tree/master?ref_type=heads

NOTE: for BGL SL_Decay has different defaults and internalparams, in particular, 
    - Different resonance masses for Blaschke factors,
      Additionally, for g, the heaviest resonance is removed as considered to close to threshold
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
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the SL_Decay output
data_BGL = np.loadtxt("checks/checks_Hammer_BdToDstFFBGL.txt", float, skiprows=1).T
# data_BGL = np.loadtxt("checks/checks_Hammer_test.txt", float, skiprows=1).T
qsq = data_BGL[0]
# hammerFF_spectrum = {
#     "f"     : data_BGL[2],
#     "g"     : data_BGL[3],
#     "F1"    : data_BGL[4],
#     "F2"    : data_BGL[5],
#     "Blf"   : data_BGL[6],
#     "Blg"   : data_BGL[7],
#     "BlP1"  : data_BGL[8],
#     "Phif"  : data_BGL[9],
#     "Phig"  : data_BGL[10],
#     "PhiF1" : data_BGL[11],
#     "PhiF2" : data_BGL[12],
#     "z"     : data_BGL[1]
# }

# hammerFF_spectrum = {
#     "f"     : data_BGL[1],
#     "g"     : data_BGL[2],
#     "F1"    : data_BGL[3],
#     "F2"    : data_BGL[4],
# }
hammerFF_spectrum = {
    "f"     : data_BGL[1],
    "g"     : data_BGL[2],
    "Fm"    : data_BGL[3],
    "Fp"    : data_BGL[4],
}

hammerFF_defaults = {
    "a0" : 0.00038,
    "a1" : 0.026905,
    "a2" : 0.,
    "b0" : 0.00055,
    "b1" : -0.0020370,
    "b2" : 0.,
    "c1" : -0.000433,
    "c2" : 0.005353,
    "d0" : 0.007,
    "d1" : -0.036
}

# Initializing slop prediction and aligning parameters
slopFF = BdToDstFF.BGL()
slopFF.set_ff(**hammerFF_defaults)
# Getting spectrum from SLOP
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff_gfF1F2_basis")
# Add hammer factor of etaEW*Vcb
etaEWVcb = slopFF.internalparams["Vcb"]*slopFF.internalparams["etaEW"]
slopFF_spectrum = {k : slopFF_spectrum[k]/etaEWVcb for k in slopFF_spectrum}

# Adjusting hammer for equivalennt basis
Mb = slopFF.internalparams["Mb"]
Mc = slopFF.internalparams["Mc"]
Mb2 = Mb**2
Mb3 = Mb**3
rC = Mc/Mb
rC2 = rC**2
w = slopFF.w(qsq)
w2 = w**2
Fpf = (rC-w)/(2.*rC*Mb2*(w2 - 1))
FpF1 = 1./(2.*rC*Mb3*(w2 - 1))
Fmf = (rC+w)/(2.*rC*Mb2*(w2 - 1))
FmF1 = 1./(2.*rC*Mb3*(w2 - 1))*(rC2-1)/(1 + rC2 - 2*rC*w)
FmP1 = np.sqrt(rC)*(rC+1)/(Mb*(1 + rC2 - 2*rC*w))
F1 = (1./FpF1)*(hammerFF_spectrum["Fp"] - Fpf*hammerFF_spectrum["f"])
F2 = (1./FmP1)*((1+rC)/np.sqrt(rC))*(hammerFF_spectrum["Fm"] - Fmf*hammerFF_spectrum["f"] - FmF1*F1)
hammerFF_spectrum["F1"] = F1
hammerFF_spectrum["F2"] = F2
hammerFF_spectrum["g"] = hammerFF_spectrum["g"]*2.

# Making comparison plot
chk.make_comparison_plot(slopFF_spectrum, hammerFF_spectrum, qsq, "Hammer v1.2.1",
    ["f", "g", "F1", "F2"],
    [r"$f$", r"$g$", r"$\mathcal{F}_1$", r"$\mathcal{F}_2$"],
    "BdToDstFFBGL", "checks/check_hammer_{}_{}.png")

# slopFF_bl = chk.get_additional_spectrum_BGL(qsq, slopFF)
# chk.make_comparison_plot(slopFF_bl, hammerFF_spectrum, qsq, "Hammer v1.2.1",
#     ["Blg", "Blf", "BlP1", "Phig", "Phif", "PhiF1", "PhiF2", "z"],
#     [],
#     "BdToDstFFBGL", "checks/check_hammer_{}_{}.png")