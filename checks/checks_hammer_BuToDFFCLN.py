"""
Cross-check with SL_Decay, see
- https://cds.cern.ch/record/2313977/files/LHCb-INT-2018-015.pdf
- https://gitlab.cern.ch/scali/SL_Decay/-/tree/master?ref_type=heads

NOTE: for this comparison, convert HAMMER basis to h basis and get SLOP predictions in h basis
"""
from slophep.Predictions.FormFactorsBToP import BuToDFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data = np.loadtxt("checks/checks_hammer_BuToDFFCLN.txt", float, skiprows=1).T
qsq = data[0]
hammerFF_spectrum = {
    "Fs"    : data[1],
    "f0"    : data[2],
    "f+"    : data[3],
}

# Initializing slop prediction
slopFF = BuToDFF.CLN()
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

ff = ["f0", "f+"]
fflabel = [r"$f_{0}$", r"$f_{+}$"]
annotation = r"""Notes:
- SM only parameterisation, $f_T = 0$
"""
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BuToDFFCLN_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.2.1")
    cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)