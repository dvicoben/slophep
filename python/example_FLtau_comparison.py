"""
Example script to see how to plot (binned) experimental results with predictions
Numbers for FL taken for https://lhcbproject.web.cern.ch/Publications/LHCbProjectPublic/LHCb-PAPER-2023-020.html

This produces three plots:
- Central values for all FF schemes with an experimental result (no errorbands)
- Binned prediction (HPQCD), with error, and experimental result
- Unbinned prediction (HPQCD) i.e. full q2 range, with error, and experimental result
"""
from slophep.Predictions.SamplingFluctuate import SamplingHelper
from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV import BdToDstFF
from slophep.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

setPlotParams()

# The experimental result from https://lhcbproject.web.cern.ch/Publications/LHCbProjectPublic/LHCb-PAPER-2023-020.html
fl_exp = np.array([0.52, 0.34])
fl_exp_err_stat = np.array([0.07, 0.08])
fl_exp_err_syst = np.array([0.04, 0.02])
fl_exp_err = np.sqrt(fl_exp_err_stat**2 + fl_exp_err_syst**2)

# A quick comparison of several FF schemes and central values with the experimental value:
obs_cln = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.CLN)
obs_blpr = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.BLPR)
obs_bgl = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.BGL)
obs_bsz = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.BSZ)
obs_hpqcd = BdToDstEllNuPrediction("tau", "tau", BdToDstFF.HPQCD)
obs = "FL"
p = obs_cln.plot_obs_prediction(obs, label="CLN")
p = obs_blpr.plot_obs_prediction(obs, label="BLPR", plot=p)
p = obs_bgl.plot_obs_prediction(obs, label="BGL", plot=p)
p = obs_bsz.plot_obs_prediction(obs, label="BSZ", plot=p)
p = obs_hpqcd.plot_obs_prediction(obs, label="HPQCD", plot=p)
# Plot the data on too:
q2_bin_edges = np.array([obs_hpqcd.q2min, 7, obs_hpqcd.q2max])
q2_bin_c = 0.5*(q2_bin_edges[1:] + q2_bin_edges[:-1])
q2_bin_w = 0.5*(q2_bin_edges[1:] - q2_bin_edges[:-1])
p[1].errorbar(q2_bin_c, fl_exp, xerr=q2_bin_w, yerr=fl_exp_err, 
              color='k', markersize=2, fmt='o', capsize=2, label="LHCb-PAPER-2023-020")
p[1].set(ylabel=r"$F_L$")
p[1].legend()
p[0].savefig("output/FLtau_comparison.png", bbox_inches='tight', dpi=150)



# Now get the binned prediction and error - will use the HPQCD FFs for this
# The central values
fl_pred = np.array([obs_hpqcd.fl_bin(iq, jq) for iq, jq in zip(q2_bin_edges[:-1], q2_bin_edges[1:])])
print(f"HPQCD Prediction: {fl_pred}")
# Fluctuations for errorbands
Nfluct = 2000
fluct = SamplingHelper(obs_hpqcd)
fluct.set_params_from_configfile("data/BToDstFF_HPQCD_COV_arXiv230403137.json")
fluct.fluctuate(Nfluct)
# Getting errorbands for the binned prediction
# This has to carry out the integration for every fluctuation - can take a bit of time
fl_pred_err = np.array([fluct.get_error("fl_bin", [iq, jq]) for iq, jq in zip(q2_bin_edges[:-1], q2_bin_edges[1:])])
fl_pred_errlo = np.array([ifl[0] for ifl in fl_pred_err])
fl_pred_errhi = np.array([ifl[1] for ifl in fl_pred_err])
# Plotting:
fig, ax = plt.subplots(1, 1)
ax.set(xlabel=r"$q^2$ [GeV$^2$]", ylabel=r"$F_L$")
ax.stairs(fl_pred_errhi, q2_bin_edges, baseline=fl_pred_errlo,
          color="blue", alpha=0.25, fill=True, label="Binned HPQCD arXiv:2304.03137")
ax.errorbar(q2_bin_c, fl_exp, xerr=q2_bin_w, yerr=fl_exp_err, 
            color='k', markersize=2, fmt='o', capsize=2, label="LHCb-PAPER-2023-020")
ax.legend()
fig.savefig("output/FLtau_bin_comparison.png", bbox_inches='tight', dpi=150)



# For a full q2-range band instead:
# If we want a band for the full q2 range, we can use this function to make things easier:
def get_spectrum_dict(qsq, obs, attr, fluct):
    res = {"val": [], "lo": [], "hi": []}
    for iq2 in qsq:
        o = getattr(obs, attr)(iq2)
        o_err = fluct.get_error(attr, [iq2])
        
        res["val"].append(o)
        res["lo"].append(o_err[0])
        res["hi"].append(o_err[1])
    return res

q2_points = np.linspace(obs_hpqcd.q2min+1e-6, obs_hpqcd.q2max-1e-6, 100)
fl_full = get_spectrum_dict(q2_points, obs_hpqcd, "fl", fluct)
# Now plotting full band + experimental result:
fig, ax = plt.subplots(1, 1)
ax.set(xlabel=r"$q^2$ [GeV$^2$]", ylabel=r"$F_L$")
ax.plot(q2_points, fl_full["val"], 'b-')
ax.fill_between(q2_points, fl_full["lo"], fl_full["hi"], 
                color='b', alpha=0.25, label="HPQCD arXiv:2304.03137",
                edgecolor="none")
ax.errorbar(q2_bin_c, fl_exp, xerr=q2_bin_w, yerr=fl_exp_err, 
            color='k', markersize=2, fmt='o', capsize=2, label="LHCb-PAPER-2023-020")
ax.legend()
fig.savefig("output/FLtau_full_comparison.png", bbox_inches='tight', dpi=150)