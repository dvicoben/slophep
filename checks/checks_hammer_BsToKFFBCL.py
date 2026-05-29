"""
Cross-check with HAMMER
"""
from slophep.Predictions.FormFactorsBToP import BsToKFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data_BGL = np.loadtxt("checks/checks_hammer_BsToKFFBCL.txt", float, skiprows=1).T
qsq = data_BGL[0]
hammerFF_spectrum = {
    "f0"     : data_BGL[1],
    "f+"     : data_BGL[2],
}


# SLOP predictions
slopFF = BsToKFF.BCL()
q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)


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