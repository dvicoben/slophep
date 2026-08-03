from slophep.Observables import BdToDstEllNuPrediction
from slophep.Tools import ErrorSampler
from slophep.FormFactors import BdToDstFF

from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

# Number of fluctuations we will use for errorbands
Nfluct = 5000

# Initialise observables we will fluctuate
obs_hpqcd_mu = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.HPQCD2023())
obs_hpqcd_tau = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.HPQCD2023())

# Both predictions use the same parameter, so only need one error sampler
fconfig_hpqcd = "data/BToDstFF_HPQCD_COV_arXiv230403137.json"
errsampler = ErrorSampler.create_from_configfile(fconfig_hpqcd)
errsampler.fluctuate(Nfluct)

# We are interested in the full q2 range, so some helper functions for that:
def get_spectrum_float(qsq: list[float], 
                       obs: BdToDstEllNuPrediction, 
                       attr: str, 
                       fluct: ErrorSampler) -> dict[str, list]:
    res = {"val" : [], "lo" : [], "hi" : []}
    for iq2 in qsq:
        o_err = fluct.get_error(obs, attr, [iq2])
        o = fluct.get_central(obs, attr, [iq2])
        res["val"].append(o)
        res["lo"].append(o_err[0])
        res["hi"].append(o_err[1])
    return res

def plot_spectrum(qsq: list[float], 
                  res_list: list[dict], 
                  label_list: list[str], 
                  obs_label: list[str], 
                  savepath: str) -> None:
    for ires, label in zip(res_list, label_list):
        lines = plt.plot(qsq, ires["val"])
        fill = plt.fill_between(qsq, ires["lo"], ires["hi"], color=lines[0].get_color(), alpha=0.15, label=label)
    
    plt.xlabel(r"$q^2$")
    plt.ylabel(obs_label)
    plt.legend()
    plt.savefig(savepath, bbox_inches='tight')
    plt.clf()
    plt.close()


# Now lets plot comparisons in the full q2 range
npoints = 100
qsq = np.linspace(obs_hpqcd_mu.q2min+1e-6, obs_hpqcd_mu.q2max-1e-6, npoints)

setPlotParams()
mu_res = get_spectrum_float(qsq, obs_hpqcd_mu, "dBRdq2", errsampler)
tau_res = get_spectrum_float(qsq, obs_hpqcd_tau, "dBRdq2", errsampler)
plot_spectrum(qsq, [mu_res, tau_res], 
              [r"$B^0 \to D^* \mu\nu$", r"$B^0 \to D^* \tau \nu$"],
              r"$\mathrm{d}\mathcal{B}/\mathrm{d}q^2$",
              "output/plot_q2spectrum_dBR_mutau.png")