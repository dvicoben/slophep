"""
Cross-check with HAMMER

NOTE: for BGL in HAMMER has different defaults and differs in: 
    - Different ffpar defaults
    - Additional scaling of 1./etaEW*Vcb
    - Different FF basis
For this comparison: 
    - We set the ffpar to hammer defaults in SLOP
    - Scale SLOP predictions according to the factor of 1./etaEW*Vcb.
    - Get SLOP predictions as f, g, F1, F2
    - Convert HAMMER basis to f, g, F1, F2
"""
from slophep.Predictions.FormFactorsBToV import BsToDsstFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data_BGL = np.loadtxt("checks/checks_hammer_BsToDsstFFBGL.txt", float, skiprows=1).T
qsq = data_BGL[0]
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

# Getting SLOP default predictions
slopFF = BsToDsstFF.BGL()
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff_gfF1F2_basis")
# Add hammer factor of etaEW*Vcb
etaEWVcb = slopFF.internalparams["Vcb"]*slopFF.internalparams["etaEW"]

slopFF_aligned = BsToDsstFF.BGL()
slopFF_aligned.set_ff(**hammerFF_defaults)
slopFF_aligned_spectrum = chk.get_spectrum_slop(slopFF_aligned, qsq, "get_ff_gfF1F2_basis")
slopFF_aligned_spectrum = {k : slopFF_aligned_spectrum[k]/etaEWVcb for k in slopFF_aligned_spectrum}

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


# Comparison plots
ff = ["f", "g", "F1", "F2"]
fflabel = [r"$f$", r"$g$", r"$\mathcal{F}_1$", r"$\mathcal{F}_2$"]
annotation = r"""Notes:
- HAMMER divides its FFs by
  a factor of $\eta_{EW} V_{cb}$
- SLOP BGL coefficients are 
  HAMMER's divided by $\eta_{EW} V_{cb}$
- For `aligned', we set BGL 
  coefficients to HAMMER defaults 
  and scale by $\eta_{EW} V_{cb}$ 
  as in HAMMER - the effect is 
  the same (lines overlap)
"""
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BsToDsstFFBGL_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_slop_prediction(qsq, slopFF_aligned_spectrum[iff], "SLOP (aligned)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.2.1")
    cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)