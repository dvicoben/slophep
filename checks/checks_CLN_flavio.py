from slophep.Predictions.FormFactorsBToV import BdToDstFF

import numpy as np
import matplotlib.pyplot as plt
import flavio

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


# Getting the SLOP predictions for CLN:
def get_spectrum(ff):
    q2max = (ff.internalparams["Mb"] - ff.internalparams["Mc"])**2
    qsq = np.linspace(0.012, q2max-1e-3, 100)
    
    res = {}
    for iq2 in qsq:
        iff = ff.get_ff(iq2)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return qsq, {k : np.array(res[k]) for k in res}


btodst_cln = BdToDstFF.CLN2()
setffpar = perpare_dst_params_from_flavio(par)
btodst_cln.set_ff(**setffpar)
btodst_qsq, btodst_clnff = get_spectrum(btodst_cln)



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

btodst_clnflavio = get_flavio_spectrum(cln_btov_flavio, "B->D*", par, btodst_qsq)


# Plotting them together
def make_comparison_plot(sloppred, otherpred, qsq, prefix):
    for ipred in sloppred:
        fig, ax = plt.subplots(1, 1)
        ax.plot(qsq, sloppred[ipred], 'b-', label="SLOP CLN")
        ax.plot(qsq, otherpred[ipred], 'r--', label="flavio CLN")
        ax.set(xlabel = r"$q^2$", ylabel=ipred, title=prefix)
        ax.legend()
        plt.savefig(f"checks/checks_flavio_{prefix}_{ipred}.png", 
                    bbox_inches = 'tight',
                    dpi=100)

make_comparison_plot(btodst_clnff, btodst_clnflavio, btodst_qsq, "BdToDstFFCLN2")