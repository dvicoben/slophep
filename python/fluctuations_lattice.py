from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.SamplingFluctuate import SamplingHelper
import slophep.Predictions.FormFactorsBToV.BdToDstFF as BdToDstFF

from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

# Number of fluctuations we will use for errorbands
Nfluct = 5000
obs_hpqcd = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.HPQCD)

# Lets do the same for the form factors:
ff_hpqcd = BdToDstFF.HPQCD()
ff_flab = BdToDstFF.BGL_FLabMILC()
ff_jlqcd = BdToDstFF.BGL_JLQCD()

ff_hpqcd_fluct = SamplingHelper(ff_hpqcd)
fconfig_hpqcd = "data/BToDstFF_HPQCD_COV_arXiv230403137.json"
ff_hpqcd_fluct.set_params_from_configfile(fconfig_hpqcd)
ff_hpqcd_fluct.fluctuate(Nfluct)
fferr_hpqcd = ff_hpqcd_fluct.get_error("get_ff", [5.0])

ff_flab_fluct = SamplingHelper(ff_flab)
fconfig_flab = "data/BToDstFF_BGL_FlabMILC_COV_arXiv210514019.json"
ff_flab_fluct.set_params_from_configfile(fconfig_flab)
ff_flab_fluct.fluctuate(Nfluct)
fferr_flab = ff_flab_fluct.get_error("get_ff", [5.0])

ff_jlqcd_fluct = SamplingHelper(ff_jlqcd)
fconfig_jlqcd = "data/BToDstFF_BGL_JLQCD_COV_arXiv230605657.json"
ff_jlqcd_fluct.set_params_from_configfile(fconfig_jlqcd)
ff_jlqcd_fluct.fluctuate(Nfluct)
fferr_jlqcd = ff_jlqcd_fluct.get_error("get_ff", [5.0])

print(ff_hpqcd.get_ff(5.0))
print(fferr_hpqcd)
print(ff_flab.get_ff(5.0))
print(fferr_flab)
print(ff_jlqcd.get_ff(5.0))
print(fferr_jlqcd)

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
        plt.savefig(f"output/plot_lattice_comp_{iobs}.png", bbox_inches='tight')
        plt.clf()
        plt.close()


# Now lets plot comparisons in the full q2 range
npoints = 100
qsq = np.linspace(obs_hpqcd.q2min+1e-6, obs_hpqcd.q2max-1e-6, npoints)

setPlotParams()

hpqcd_ff_res = get_spectrum_dict(qsq, ff_hpqcd, "get_ff_gfF1F2_basis", ff_hpqcd_fluct)
flab_ff_res = get_spectrum_dict(qsq, ff_flab, "get_ff_gfF1F2_basis", ff_flab_fluct)
jlqcd_ff_res = get_spectrum_dict(qsq, ff_jlqcd, "get_ff_gfF1F2_basis", ff_jlqcd_fluct)
plot_spectrum_dict(qsq, ["g", "f", "F1", "F2"], 
                   [hpqcd_ff_res, flab_ff_res, jlqcd_ff_res], 
                   ["HPQCD arXiv:2304.03137", "Flab-MILC arXiv:2105.14019", "JLQCD arXiv:2306.05657"],
                   [r"$g$", r"$f$", r"$\mathcal{F}_1$", r"$\mathcal{F}_2$"])