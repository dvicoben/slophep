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


slopFF = BdToDstFF.CLN2()
setffpar = perpare_dst_params_from_flavio(par)
slopFF.set_ff(**setffpar)

q2max = (slopFF.internalparams["Mb"] - slopFF.internalparams["Mc"])**2
q2min = slopFF.par["m_mu"]**2
qsq = np.linspace(q2min+1e-6, q2max-1e-6, 100)

slopFF_spectrum = chk.get_spectrum_slop(slopFF, qsq, "get_ff")

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

flavioFF_specturm = get_flavio_spectrum(cln_btov_flavio, "B->D*", par, qsq)


# Plotting them together
chk.make_comparison_plot(slopFF_spectrum, flavioFF_specturm, qsq, "flavio",
    ["A0", "A1", "A12", "V", "T1", "T2", "T23"],
    [r"$A_0$", r"$A_1$", r"$A_{12}$", r"$V$", r"$T_1$", r"$T_2$", r"$T_{23}$"],
    "BdToDstFFCLN2", "checks/check_flavio_{}_{}.png")