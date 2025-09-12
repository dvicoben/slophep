from slophep.Predictions.FormFactorsBToV import BdToDstFF
from slophep.utils import setPlotParams
import check_utils as chk
import numpy as np
import matplotlib.pyplot as plt
import flavio

setPlotParams()

par = flavio.default_parameters.get_central_all()
def perpare_dst_params_from_flavio(par):
    """A small helper function so we can make sure to set the 
       same parameters in SLOP as what flavio is using
    """
    ffpar = {
        "RhoSq" : par["B->D* CLN rho2"],
        "h_A1"  : par["B->D* CLN h_A1(1)"],
        "R1"    : par["B->D* CLN R_1(1)"],
        "R2"    : par["B->D* CLN R_2(1)"],
    }
    mB, mV = par["m_B0"], par["m_D*+"]
    w0 = (mB**2 + mV**2) / (2*mB*mV)
    R2_0 = ffpar["R2"] + 0.11*(w0-1) - 0.06*(w0-1)**2 # R_2 at q^2=0
    R0_0 = (mB + mV - mB * R2_0 + mV * R2_0)/(2 * mV) # R_0 at q^2=0
    R0 = R0_0 - (- 0.11*(w0-1) + 0.01*(w0-1)**2)
    ffpar["R0"] = R0
    return ffpar

# Using the SLOP defaults
slopFF = BdToDstFF.CLN()

# Aligning parameters
slopFF_aligned = BdToDstFF.CLN2()
setffpar = perpare_dst_params_from_flavio(par)
slopFF_aligned.set_ff(**setffpar)

q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")
slopFF_aligned_spectrum = chk.get_spectrum_slop(slopFF_aligned, qsq, "get_ff")

# Getting the flavio predictions:
import flavio
from flavio.physics.bdecays.formfactors.b_v.clnexp import ff as cln_btov_flavio
# For some reason, in https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py
# flavio sets B->D* as B0->D*0 rather than B0->D*+ but then the isgur_wise function used for 
# the tensor ffs uses the correct B0->D*+
# We need this small workaround to ensure things are being done correctly!
process_dict = {
    "B->D*" : {'B': 'B0', 'V': 'D*+',   'q': 'b->c'}
}
flavio.physics.bdecays.formfactors.b_v.clnexp.process_dict = process_dict

def get_flavio_spectrum(fffunc, process, par, qsq):
    res = {}
    for iq2 in qsq:
        iff = fffunc(process, iq2, par, 4.8)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return {k : np.array(res[k]) for k in res}

flavioFF_spectrum = get_flavio_spectrum(cln_btov_flavio, "B->D*", par, qsq)

# Comparison plots
ff = ["A0", "A1", "A12", "V", "T1", "T2", "T23"]
fflabel = [r"$A_0$", r"$A_1$", r"$A_{12}$", r"$V$", r"$T_1$", r"$T_2$", r"$T_{23}$"]
annotation = r"""Notes:
- For `aligned', we set $\rho^2$, $h_{A1}(1)$, 
  $R_1(1)$, $R_2(1)$, and $R_0(1)$ to flavio 
  defaults
"""
extra_note = "- SLOP CLN is SM only so\n" +r"  $T_1 = T_2 = T_{23} = 0$, flavio behaviour"+"\n"+r"  can be reproduced with CLN2"
for iff, ifflabel in zip(ff, fflabel):
    savepath = f"checks/check_flavio_BdToDstFFCLN_{iff}.png"
    cplot = chk.ComparisonPlot(ifflabel)
    cplot.add_slop_prediction(qsq, slopFF_spectrum[iff], "SLOP (default)")
    cplot.add_slop_prediction(qsq, slopFF_aligned_spectrum[iff], "SLOP (aligned)")
    cplot.add_comparison_prediction(qsq, flavioFF_spectrum[iff], "flavio")
    note = annotation
    if "T" in iff:
        note += extra_note
    cplot.annotate(note, 1.01, 0.5)
    cplot.makeplot()
    cplot.savefig(savepath)