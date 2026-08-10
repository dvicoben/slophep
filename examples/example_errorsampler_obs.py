from slophep.Observables import BdToDstEllNuPrediction
from slophep.Tools import ErrorSampler
from slophep.FormFactors import BdToDstFF

from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

# Number of fluctuations we will use for errorbands
Nfluct = 2000

# Initialise observables we will fluctuate
obs_hpqcd = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.HPQCD2023())
obs_bsz = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.BSZ())

# Generate parameter fluctations
hpqcd_errfluct = ErrorSampler.create_from_configfile("data/BdToDstFF_HPQCD_COV_arXiv230403137.json")
hpqcd_errfluct.fluctuate(Nfluct)
bsz_errfluct = ErrorSampler.create_from_configfile("data/BdToDstFF_BSZ_COV_arXiv181100983.json")
bsz_errfluct.fluctuate(Nfluct)

# Print nominal values and errors at particular q2 value
print(hpqcd_errfluct.get_central(obs_hpqcd, "J", [5.0]))
print(hpqcd_errfluct.get_error(obs_hpqcd, "J", [5.0]))
print(bsz_errfluct.get_central(obs_bsz, "J", [5.0]))
print(bsz_errfluct.get_error(obs_bsz, "J", [5.0]))

# We are interested in the full q2 range, so some helper functions for that:
def get_spectrum_dict(qsq: list[float], 
                      obs, 
                      attr: str,
                      fluct: ErrorSampler):
    res = {}
    for iq2 in qsq:
        o_err = fluct.get_error(obs, attr, [iq2])
        o = fluct.get_central(obs, attr, [iq2])
        for iobs in o:
            if iobs not in res:
                res[iobs] = {"val" : [], "lo" : [], "hi" : []}
            res[iobs]["val"].append(o[iobs])
            res[iobs]["lo"].append(o_err[iobs][0])
            res[iobs]["hi"].append(o_err[iobs][1])
    return res

def plot_spectrum_dict(qsq: list[float], 
                       obslist: list[str], 
                       res_list: list[dict], 
                       label_list: list[str], 
                       obs_labels: list[str],):
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
        plt.savefig(f"output/example_errorsampler_obs_{iobs}.png", bbox_inches='tight')
        plt.clf()
        plt.close()


# Now lets plot comparisons in the full q2 range
npoints = 100
qsq = np.linspace(obs_hpqcd.q2min+1e-6, obs_hpqcd.q2max-1e-6, npoints)

setPlotParams()

hpqcd_J_res = get_spectrum_dict(qsq, obs_hpqcd, "J", hpqcd_errfluct)
bsz_J_res = get_spectrum_dict(qsq, obs_bsz, "J", bsz_errfluct)
plot_spectrum_dict(qsq, [ielem for ielem in hpqcd_J_res], 
                   [hpqcd_J_res, bsz_J_res], 
                   ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"],
                   [r"$J_{"+str(ielem)+r"}$" for ielem in hpqcd_J_res])

hpqcd_Juni_res = get_spectrum_dict(qsq, obs_hpqcd, "uniang_obs", hpqcd_errfluct)
bsz_Juni_res = get_spectrum_dict(qsq, obs_bsz, "uniang_obs", bsz_errfluct)
plot_spectrum_dict(qsq, ["FL", "AFB", "FLt"], 
                   [hpqcd_Juni_res, bsz_Juni_res], 
                   ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"],
                   [r"$F_L$", r"$A_{FB}$", r"$\tilde{F}_L$"])