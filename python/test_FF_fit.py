"""
Simple script to test FF fits using BToDstObs prediction as fit model

NOTE: At the moment fits are quite slow as:
    - Method (as implemented) recomputes the angular integrals at every step 
      (can be reduced to compute only once, since these are always the same)
    - Need to recompute and integrate J_i at every fit step (cannot really be remedied)
"""

import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uncert
from iminuit import Minuit

from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.FormFactorsBToDst import BGL

obs = BToDstEllNuPrediction("mu", "mu", BGL)
# { # https://arxiv.org/pdf/1707.09509 (Third column in table V)
#     "a0" : 0.0209,
#     "a1" : 0.33,
#     "a2" : 0.6,
#     "b0" : 0.01218,
#     "b1" : 0.046,
#     "b2" : 0.48,
#     "c1" : 0.0063,
#     "c2" : 0.062,
#     "d0" : 0.0595,
#     "d1" : -0.218
# }
ff_gen = { # https://arxiv.org/pdf/1707.09509 (Third column in table V)
    "a0" : 0.0209,
    "a1" : 0.33,
    "a2" : 0.0,
    "b0" : 0.01218,
    "b1" : 0.046,
    "b2" : 0.0,
    "c1" : 0.0063,
    "c2" : 0.0,
    "d0" : 0.0595,
    "d1" : 0.0
}
obs.set_ff(ff_gen)
print("Generated Values:")
print(obs.FF.ffpar)

ctd_bins = 4
ctl_bins = 4
chi_bins = 5
# Gets a histogram of the PDF with specified binning
h, b, angint, j_bins = obs.PDF_hist(10, ctd_bins, ctl_bins, chi_bins)

N = 2e6
h_norm = h/np.sum(h)
h_data = N*h_norm
# Expect non-poisson errors + unfolding will increase uncertainty - set to 2*sigma_Poisson
h_data_err = 2.0*np.sqrt(h_data)


pred = BToDstEllNuPrediction("mu", "mu", BGL)


def calc_pdf(ff: list[float]):
    # ffd = {
    #     "a0" : ff[0],
    #     "a1" : ff[1],
    #     "a2" : ff[2],
    #     "b0" : 0.01218,
    #     "b1" : ff[3],
    #     "b2" : ff[4],
    #     "c1" : ff[5],
    #     "c2" : ff[6],
    #     "d0" : ff[7],
    #     "d1" : ff[8],
    # }
    ffd = {
        "a0" : ff[0],
        "a1" : ff[1],
        "a2" : 0.0,
        "b0" : 0.01218,
        "b1" : ff[2],
        "b2" : 0.0,
        "c1" : ff[3],
        "c2" : 0.0,
        "d0" : ff[4],
        "d1" : 0.0,
    }
    pred.set_ff(ffd)
    print(pred.FF.ffpar)
    hist, *_ = pred.PDF_hist(*b)
    return hist

def fom(ff: list[float]):
    # Simple chi2 fit - assume uncorrletad bin yields but this will likely not be the case
    pred = calc_pdf(ff)
    h_pred = np.sum(h_data)*(pred/np.sum(pred))
    chiSq = np.sum((h_data - h_pred)**2/(h_data_err**2))
    # invcov = np.linalg.inv(cov)
    # delta_p = (h_data - h_pred) # may need to unroll/flatten histograms
    # chiSq = np.dot(delta_p, np.dot(invcov, delta_p))
    return chiSq


# fitpars = ["a0", "a1", "a2", "b1", "b2", "c1", "c2", "d0", "d1"]
fitpars = ["a0", "a1", "b1", "c1", "d0"]
init_vals = [obs.FF.ffpar[ipar] for ipar in fitpars]
m = Minuit(fom, init_vals)
m.migrad()
m.hesse()

print("Function Minimum:", m.fval)
if m.valid:
    print("Fit converged")
else:
    print("Fit did not converge")

if m.accurate:
    print("Covariance matrix accurate")
else:
    print("Covariance matrix not accurate")

reslist = [elem.value for elem in m.params]
errlist = [elem.error for elem in m.params]
print(fitpars)
print(init_vals)
print(reslist)
print(errlist)

import json
res = {
    "val" : reslist,
    "err" : errlist,
    "pars" : fitpars,
    "gen" : ff_gen,
}
outpath = f"output/fitres_FFBGL1110_N2e6_B{ctd_bins}_{ctl_bins}_{chi_bins}.json"
with open(outpath, "w") as outfile:
    json.dump(res, outfile, indent=2)