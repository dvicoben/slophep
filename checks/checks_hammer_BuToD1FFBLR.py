"""
Cross-check with SL_Decay, see
- https://cds.cern.ch/record/2313977/files/LHCb-INT-2018-015.pdf
- https://gitlab.cern.ch/scali/SL_Decay/-/tree/master?ref_type=heads

NOTE: for this comparison, convert HAMMER basis to h basis and get SLOP predictions in h basis
"""
from slophep.Predictions.FormFactorsBToDstst import BuToD1FF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np

setPlotParams()

# load in the Hammer output
data = np.loadtxt("checks/checks_hammer_BuToD1FFBLR.txt", float, skiprows=1).T
qsq = data[0]
hammerFF_spectrum = {
    "fS"    : data[1],
    "fV1"   : data[2],
    "fV2"   : data[3],
    "fV3"   : data[4],
    "fA"    : data[5],
    "fT1"   : data[6],
    "fT2"   : data[7],
    "fT3"   : data[8]
}

# Initializing slop prediction and aligning parameters
slopFF = BuToD1FF.BLR()
# Getting spectrum from SLOP
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

ff = ["fV1", "fV2", "fV3", "fA", "fT1", "fT2", "fT3"]
fflabel = [r"$f_{V_1}$", r"$f_{V_2}$", r"$f_{V_3}$", r"$f_{A}$", r"$f_{T_1}$", r"$f_{T_2}$", r"$f_{T_3}$"]
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BuToD1FFBLR_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.2.1")
    cplot.makeplot()
    cplot.savefig(savepath)