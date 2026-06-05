"""
Cross-check with HAMMER
"""
from slophep.Predictions.FormFactorsBToP import BsToKFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data = np.loadtxt("checks/checks_hammer_BsToKFFBCL.txt", float, skiprows=1).T
qsq = data[0]
hammerFF_spectrum = {
    "f0"     : data[1],
    "f+"     : data[2],
}


# SLOP predictions
slopFF = BsToKFF.BCL()
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")


ff = ["f+", "f0"]
fflabel = [r"$f_+$", r"$f_0$"]
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BsToKFFBCL_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default, BCL)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.4.1")
    # cplot.annotate(annotation, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)