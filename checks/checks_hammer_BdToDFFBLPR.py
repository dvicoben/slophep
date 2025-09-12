"""
Cross-check with SL_Decay, see
- https://cds.cern.ch/record/2313977/files/LHCb-INT-2018-015.pdf
- https://gitlab.cern.ch/scali/SL_Decay/-/tree/master?ref_type=heads

NOTE: for this comparison, convert HAMMER basis to h basis and get SLOP predictions in h basis
"""
from slophep.Predictions.FormFactorsBToP import BdToDFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data = np.loadtxt("checks/checks_hammer_BdToDFFBLPR.txt", float, skiprows=1).T
qsq = data[0]
hammerFF_spectrum = {
    "Fs"    : data[1],
    "f0"    : data[2],
    "f+"    : data[3],
    "fT"    : data[4],
}

# Initializing slop prediction
slopFF = BdToDFF.BLPR()
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

# Initializing slop prediction and aligning parameters
slopFF_aligned = BdToDFF.BLPR()
slopFF_aligned_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")
slopFF_aligned_spectrum["fT"] = slopFF_aligned_spectrum["fT"]/(slopFF_aligned.internalparams["Mc"] + slopFF_aligned.internalparams["Mb"])

ff = ["f0", "f+", "fT"]
fflabel = [r"$f_{0}$", r"$f_{+}$", r"$f_{T}$"]
annotation = r"""Notes:
- SLOP follows flavio conventions
  for the FFs as flavio is used to
  calculate observables. This differs
  only in that 
      $f_T = \frac{(m_D m_B)h_T}{2\sqrt{m_D m_B}}$
  (i.e. eqn. B4 in arXiv:1908.09398),
  rather than HAMMER's
      $f_T = \frac{h_T}{2\sqrt{m_D m_B}}$
- The `aligned' simply accounts for
  the factor of $(m_D + m_B)$
"""
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BdToDFFBLPR_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_slop_prediction(qsq, slopFF_aligned_spectrum[iff], "SLOP (aligned)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.2.1")
    cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)