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
from bd2dstlnu.Predictions.BToDstFFBLPR import BLPR


obs = BToDstEllNuPrediction("mu", "mu", BLPR)
# Gets a histogram of the PDF with specified binning
h, b, angint, j_bins = obs.PDF_hist(10, 4, 4, 5)


N = 1e6
h_norm = h/np.sum(h)
h_data = N*h_norm
# Expect non-poisson errors + unfolding will increase uncertainty - set to 2*sigma_Poisson
h_data_err = 2.0*np.sqrt(h_data)


pred = BToDstEllNuPrediction("mu", "mu", BLPR)
# {
#     "RhoSq" : 1.24,
#     "Chi21" : -0.06,
#     "Chi2p" : 0.0,
#     "Chi3p" : 0.05,
#     "Eta1"  : 0.30,
#     "Etap"  : -0.05,
#     "dV20"  : 0.0
# }

def calc_pdf(ff: list[float]):
    # This roundabout input is to fix Chi3p and dV20 for which there is no sensitivity
    ffpars = np.concatenate((ff[:3], [0.05], ff[3:], [0.0]))
    pred.set_ff_fromlist(ffpars)
    print(pred.FF.ffpar)
    hist, *_ = pred.PDF_hist(*b)
    x = np.sum(np.where(hist < 0, 1, 0))
    if x > 0:
        print(f"{x} -ve bins from {ff}")
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


fitpars = ["RhoSq", "Chi21", "Chi2p", "Eta1", "Etap"]
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