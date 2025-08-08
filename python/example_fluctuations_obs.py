from slophep.Predictions.Observables import BToDstEllNuPrediction
from slophep.Predictions.SamplingFluctuate import SamplingHelper
from slophep.Predictions.FormFactorsBToV.BToDstFF import HPQCD, BSZ

from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

# Number of fluctuations we will use for errorbands
Nfluct = 5000

wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.0,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.0
}
# Initialise observables we will fluctuate
obs_hpqcd = BToDstEllNuPrediction("mu", "mu", HPQCD)
obs_hpqcd.set_wc(wcoeffs)
obs_bsz = BToDstEllNuPrediction("mu", "mu", BSZ)
obs_bsz.set_wc(wcoeffs)

# Fluctuate at particular value of q2
obs_hpqcd_fluct = SamplingHelper(obs_hpqcd)
fconfig_hpqcd = "data/BToDstFF_HPQCD_COV_arXiv230403137.json"
obs_hpqcd_fluct.set_params_from_configfile(fconfig_hpqcd)
obs_hpqcd_fluct.fluctuate(Nfluct)
jerr_hpqcd = obs_hpqcd_fluct.get_error("J", [5.0])

obs_bsz_fluct = SamplingHelper(obs_bsz)
fconfig_bsz = "data/BToDstFF_BSZ_COV_arXiv181100983.json"
obs_bsz_fluct.set_params_from_configfile(fconfig_bsz)
obs_bsz_fluct.fluctuate(Nfluct)
jerr_bsz = obs_bsz_fluct.get_error("J", [5.0])

# Print nominal values and errors
print(obs_hpqcd.J(5.0))
print(jerr_hpqcd)
print(obs_bsz.J(5.0))
print(jerr_bsz)



# Lets do the same for the form factors:
ff_hpqcd = HPQCD()
ff_bsz = BSZ()

ff_hpqcd_fluct = SamplingHelper(ff_hpqcd)
ff_hpqcd_fluct.set_params_from_configfile(fconfig_hpqcd)
ff_hpqcd_fluct.fluctuate(Nfluct)
fferr_hpqcd = ff_hpqcd_fluct.get_error("get_ff", [5.0])


ff_bsz_fluct = SamplingHelper(ff_bsz)
ff_bsz_fluct.set_params_from_configfile(fconfig_bsz)
ff_bsz_fluct.fluctuate(Nfluct)
fferr_bsz = ff_bsz_fluct.get_error("get_ff", [5.0])
print(ff_hpqcd.get_ff(5.0))
print(fferr_hpqcd)
print(ff_bsz.get_ff(5.0))
print(fferr_bsz)


# We are interested in the full q2 range, so some helper functions for that:
def get_spectrum_dict(qsq, obs, attr, fluct):
    res = {}
    for iq2 in qsq:
        o = getattr(obs, attr)(iq2)
        o_err = fluct.get_error(attr, [iq2])
        for iobs in o:
            if iobs not in res:
                res[iobs] = {"val" : [], "lo" : [], "hi" : []}
            res[iobs]["val"].append(o[iobs])
            res[iobs]["lo"].append(o_err[iobs][0])
            res[iobs]["hi"].append(o_err[iobs][1])
    return res

def plot_spectrum_dict(qsq, obslist, res_list, label_list, obs_labels):
    
    for iobs, ilobs in zip(obslist, obs_labels):
        for res, label in zip(res_list, label_list):
            if iobs not in res:
                continue
            ires = res[iobs]
            lines = plt.plot(qsq, ires["val"])
            fill = plt.fill_between(qsq, ires["lo"], ires["hi"], color=lines[0].get_color(), alpha=0.15, label=label)
        
        plt.xlabel(r"$q^2$")
        plt.ylabel(ilobs)
        plt.legend()
        plt.savefig(f"output/plot_q2sepctrum_{iobs}.png", bbox_inches='tight')
        plt.clf()
        plt.close()


# Now lets plot comparisons in the full q2 range
npoints = 100
qsq = np.linspace(obs_hpqcd.q2min+1e-6, obs_hpqcd.q2max-1e-6, npoints)

setPlotParams()

hpqcd_ff_res = get_spectrum_dict(qsq, ff_hpqcd, "get_ff", ff_hpqcd_fluct)
bsz_ff_res = get_spectrum_dict(qsq, ff_bsz, "get_ff", ff_bsz_fluct)
plot_spectrum_dict(qsq, ["A0", "A1", "A12", "V", "T1", "T2", "T23"], 
                   [hpqcd_ff_res, bsz_ff_res], 
                   ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"],
                   [r"$A_0$", r"$A_1$", r"$A_{12}$", r"$V$", r"$T_1$", r"$T_2$", r"$T_{23}$"])

hpqcd_J_res = get_spectrum_dict(qsq, obs_hpqcd, "J", obs_hpqcd_fluct)
bsz_J_res = get_spectrum_dict(qsq, obs_bsz, "J", obs_bsz_fluct)
plot_spectrum_dict(qsq, [ielem for ielem in hpqcd_J_res], 
                   [hpqcd_J_res, bsz_J_res], 
                   ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"],
                   [r"$J_{"+str(ielem)+r"}$" for ielem in hpqcd_J_res])

hpqcd_Juni_res = get_spectrum_dict(qsq, obs_hpqcd, "uniang_obs", obs_hpqcd_fluct)
bsz_Juni_res = get_spectrum_dict(qsq, obs_bsz, "uniang_obs", obs_bsz_fluct)
plot_spectrum_dict(qsq, ["FL", "AFB", "FLt"], 
                   [hpqcd_Juni_res, bsz_Juni_res], 
                   ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"],
                   [r"$F_L$", r"$A_{FB}$", r"$\tilde{F}_L$"])