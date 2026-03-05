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
    "S"    : data[1],
    "V1"   : data[2],
    "V2"   : data[3],
    "V3"   : data[4],
    "A"    : data[5],
    "T1"   : data[6],
    "T2"   : data[7],
    "T3"   : data[8]
}

# Initializing slop prediction and aligning parameters
slopFF = BuToD1FF.BLR()
# Getting spectrum from SLOP
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

ff = ["V1", "V2", "V3", "A", "T1", "T2", "T3"]
# fflabel = [r"$h_{A1}$", r"$h_{A2}$", r"$h_{A3}$", r"$h_{V}$", r"$h_{T1}$", r"$h_{T2}$", r"$h_{T3}$"]
fflabel = ff
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_hammer_BuToD1FFBLR_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_comparison_prediction(qsq, hammerFF_spectrum[iff], "Hammer v1.2.1")
    cplot.makeplot()
    cplot.savefig(savepath)