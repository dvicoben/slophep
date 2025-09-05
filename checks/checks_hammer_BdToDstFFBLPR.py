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

# load in the Hammer output
data = np.loadtxt("checks/checks_Hammer_BdToDstFFBLPR.txt", float, skiprows=1).T
qsq = data[0]
hammerFF_spectrum = {
    "Fs"    : data[1],
    "Ff"    : data[2],
    "Fg"    : data[3],
    "Fm"    : data[4],
    "Fp"    : data[5],
    "Fzt"   : data[6],
    "Fmt"   : data[7],
    "Fpt"   : data[8]
}

# Initializing slop prediction and aligning parameters
slopFF = BdToDstFF.BLPR()
# Getting spectrum from SLOP
slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff_h_basis")

# Adjusting hammer for equivalent basis
def hammer_to_h_basis(Sqq, ham, slopFF):
    ff = {}
    
    Mb = slopFF.internalparams["Mb"]
    Mc = slopFF.internalparams["Mc"]
    sqMbMc = np.sqrt(Mb*Mc)
    Mb3 = Mb**3
    
    ff["hA1"] = (2. * sqMbMc)*ham["Ff"]/((Mb + Mc)*(Mb + Mc) - Sqq)
    ff["hV"] = (2. * sqMbMc)*ham["Fg"]
    ff["hA2"] = -(ham["Fm"] + ham["Fp"])/np.sqrt(Mc / Mb3)
    ff["hA3"] = 0.5*(2. * sqMbMc)*(ham["Fm"] - ham["Fp"])
    ff["hT3"] = ham["Fzt"]*(2. * pow(Mb * Mc, 1.5))/Mc
    ff["hT2"] = (2. * sqMbMc)*((Mb+Mc)*ham["Fmt"] + (Mb-Mc)*ham["Fpt"])/((Mb-Mc)**2 - (Mb+Mc)**2)
    ff["hT1"] = (ham["Fmt"] + (ff["hT2"] * (Mb + Mc)) / (2. * sqMbMc))*((2. * sqMbMc)/(Mb - Mc))
    return ff

hammerFF_spectrum_h = hammer_to_h_basis(qsq, hammerFF_spectrum, slopFF)

# Making comparison plot
chk.make_comparison_plot(slopFF_spectrum, hammerFF_spectrum_h, qsq, "Hammer v1.2.1",
    ["hA1", "hA2", "hA3", "hV", "hT1", "hT2", "hT3"],
    [r"$h_{A1}$", r"$h_{A2}$", r"$h_{A3}$", r"$h_{V}$", r"$h_{T1}$", r"$h_{T2}$", r"$h_{T3}$"],
    "BdToDstFFBLPR", "checks/check_hammer_{}_{}.png")
