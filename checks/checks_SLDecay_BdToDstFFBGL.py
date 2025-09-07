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
data_BGL = np.loadtxt("checks/checks_SLDecay_BdToDstFFBGL.txt", float, skiprows=24).T
qsq = data_BGL[0]
SLDecayFF_spectrum = {
    "g"     : data_BGL[2],
    "f"     : data_BGL[3],
    "F1"    : data_BGL[4],
    "F2"    : data_BGL[5],
    "Blf"   : data_BGL[6],
    "Blg"   : data_BGL[7],
    "BlP1"  : data_BGL[8],
    "Phig"  : data_BGL[9],
    "Phif"  : data_BGL[10],
    "PhiF1" : data_BGL[11],
    "PhiF2" : data_BGL[12],
    "z"     : data_BGL[1]
}

# The SLDecay default parameters
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
    "Mb"         : 5.27963,
    "chim"       : 3.894e-4,
    "chip"       : 5.131e-4,
    "chimL"      : 1.9421e-2,
    "BcStatesf"  : np.array([6.739, 6.750, 7.145, 7.150]),  # GeV
    "BcStatesg"  : np.array([6.329, 6.920, 7.020]),         # GeV
    "BcStatesP1" : np.array([6.275, 6.842, 7.250])          # GeV
}

# SLOP prediction
slopFF = BdToDstFF.BGL()
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff_gfF1F2_basis")

# Aligned SLOP prediction
slopFF_aligned = BdToDstFF.BGL()
slopFF_aligned.set_ff(**SLDecay_BGL_ffpar)
slopFF_aligned._internalparams.update(SLDecay_BGL_internalparams)
rC = slopFF_aligned.internalparams["Mc"]/slopFF_aligned.internalparams["Mb"]
slopFF_aligned_spectrum = chk.get_spectrum_slop(slopFF_aligned, qsq, "get_ff_gfF1F2_basis")
# Adjusting for equivalence with SL_Decay in F2, not sure if intended
slopFF_aligned_spectrum["F2"] = slopFF_aligned_spectrum["F2"]*np.sqrt(rC)/(1+rC)


# Comparison plots
ff = ["f", "g", "F1", "F2"]
fflabel = [r"$f$", r"$g$", r"$\mathcal{F}_1$", r"$\mathcal{F}_2$"]
annotation = r"""Notes:
- SLOP BGL coefficients are from 
  arXiv:1707.09509 (Third column in Table V)
- For `aligned' do the following:
  - Set BGL coefficients to SL Decay defaults
  - Set resonances for Blaschke factors as in
    SL Decay - in particular for $B^0 \to D^{*}$
    SL Decay has one less resonance for $g$
  - Set $\chi$'s to SL Decay values
  - $\mathcal{F}_2$ is scaled by $\sqrt{m_{B^0}/m_{D^{*}}}/(1+m_{B^0}/m_{D^{*}})$
    to match the factor in SL Decay
"""
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_SLDecay_BdToDstFFBGL_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_slop_prediction(qsq, slopFF_aligned_spectrum[iff], "SLOP (aligned)")
    cplot.add_comparison_prediction(qsq, SLDecayFF_spectrum[iff], "SL Decay")
    cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)