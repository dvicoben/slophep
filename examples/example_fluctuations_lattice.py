from slophep.Observables import BdToDstEllNuPrediction
from slophep.Tools import ErrorSampler
from slophep.FormFactors import BdToDstFF

from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

# Number of fluctuations we will use for errorbands
Nfluct = 2000
obs_hpqcd = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.HPQCD2023())

# Lets do the same for the form factors:
ff_hpqcd    = BdToDstFF.HPQCD2023()
ff_fnalmilc = BdToDstFF.BGL_FNALMILC2021()
ff_jlqcd    = BdToDstFF.BGL_JLQCD2023()

hpqcd_errfluct = ErrorSampler.create_from_configfile("data/BToDstFF_HPQCD_COV_arXiv230403137.json")
hpqcd_errfluct.fluctuate(Nfluct)
fnalmilc_errfluct = ErrorSampler.create_from_configfile("data/BToDstFF_BGL_FNALMILC_COV_arXiv210514019.json")
fnalmilc_errfluct.fluctuate(Nfluct)
jlqcd_errfluct = ErrorSampler.create_from_configfile("data/BToDstFF_BGL_JLQCD_COV_arXiv230605657.json")
jlqcd_errfluct.fluctuate(Nfluct)

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
        plt.savefig(f"output/plot_lattice_comp_{iobs}.png", bbox_inches='tight')
        plt.clf()
        plt.close()


# Now lets plot comparisons in the full q2 range
npoints = 100
qsq = np.linspace(obs_hpqcd.q2min+1e-6, obs_hpqcd.q2max-1e-6, npoints)

setPlotParams()

hpqcd_ff_res = get_spectrum_dict(qsq, ff_hpqcd   , "get_ff_gfF1F2_basis", hpqcd_errfluct)
flab_ff_res  = get_spectrum_dict(qsq, ff_fnalmilc, "get_ff_gfF1F2_basis", fnalmilc_errfluct)
jlqcd_ff_res = get_spectrum_dict(qsq, ff_jlqcd   , "get_ff_gfF1F2_basis", jlqcd_errfluct)
plot_spectrum_dict(qsq, ["g", "f", "F1", "F2"], 
                   [hpqcd_ff_res, flab_ff_res, jlqcd_ff_res], 
                   ["HPQCD arXiv:2304.03137", "Flab-MILC arXiv:2105.14019", "JLQCD arXiv:2306.05657"],
                   [r"$g$", r"$f$", r"$\mathcal{F}_1$", r"$\mathcal{F}_2$"])