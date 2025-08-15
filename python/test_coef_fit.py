"""
Simple script for angular coefficient fits

NOTE: Simplified set-up here assumes no correlations between bin yields
      In practice the unfolding means there will be some non-zero correlations
      between the yields. There may be correlations between different q2 bins
      so they may need to be fit simultaneously rather than alone as is done at the moment.
"""
import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uncert
from iminuit import Minuit

from slophep.Predictions.Observables import BToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV.BToDstFF import BLPR
from slophep.Predictions.Math.BToVMathTools import calc_norm_j


def makeJfitres(res, cov, order):
    j = uncert.correlated_values(res, cov)
    I1c = (4 - 6*j[0] + j[2] + 2*j[1])/3
    J = {order[i] : j[i] for i in range(len(j))}
    J['1c'] = I1c
    return J

obs = BToDstEllNuPrediction("mu", "mu", BLPR)
ctd_bins = 4
ctl_bins = 4
chi_bins = 5
# Gets a histogram of the PDF with specified binning
h, b, angint, j_bins = obs.PDF_hist(10, ctd_bins, ctl_bins, chi_bins)

# Scale to approx expected yields
N = 1e6
h_norm = h/np.sum(h)
h_scale = N*h_norm
# Expect non-poisson errors + unfolding will increase uncertainty - set to 2*sigma_Poisson
h_scale_err = 2.0*np.sqrt(h_scale)
# h_scale_err = 0.06*h_scale



bin_i = 0
j_i = calc_norm_j(j_bins[bin_i])
h_data = h_scale[bin_i]
h_data_err = h_scale_err[bin_i]
angarr = angint["asarray"]
jorder = angint["order"]

# ['1s', '1c', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]
def calc_pdf(j: list[float]):
    I1c = (4 - 6*j[0] + j[2] + 2*j[1])/3
    j_vec = np.concatenate(([j[0], I1c], j[1:]))
    hist = np.dot(angarr, j_vec)*(9/(32*np.pi))
    x = np.sum(np.where(hist < 0, 1, 0))
    if x > 0:
        print(f"{x} -ve bins from {j}")
    return hist

def fom(j: list[float]):
    # Simple chi2 fit - assume uncorrletad bin yields but this will likely not be the case
    pred = calc_pdf(j)
    h_pred = np.sum(h_data)*(pred/np.sum(pred))
    chiSq = np.sum((h_data - h_pred)**2/(h_data_err**2))
    # invcov = np.linalg.inv(cov)
    # delta_p = (h_data - h_pred) # may need to unroll/flatten histograms
    # chiSq = np.dot(delta_p, np.dot(invcov, delta_p))
    return chiSq


fitobs = ['1s', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]
jres = []

# Iterate over bins and fit each bin separately
# Use a chi2 which assumes uncorrelated - likely not the case so will have to incorporate covariance matrix
# Unfolding may produce correlations between different q2 bins content in which case they may need to be fitted simultaneously!
for ibin in range(len(h)):
    
    bin_i = ibin
    j_i = calc_norm_j(j_bins[bin_i])
    h_data = h_scale[bin_i]
    h_data_err = h_scale_err[bin_i]
    init_vals = [j_i[k] for k in fitobs]
    m = Minuit(fom, init_vals)
    for ipar in m.params:
        m.limits[str(ipar.name)] = (-2, 2)
    m.migrad()
    m.hesse()
    # m.minos
    print("Function Minimum:", m.fval)
    if m.valid:
        print("Fit converged")
    else:
        print("Fit did not converge")

    if m.accurate:
        print("Covariance matrix accurate")
    else:
        print("Covariance matrix not accurate")

    for i in range(len(m.params)):
        print(ibin, fitobs[i], j_i[fitobs[i]], init_vals[i], m.params[i].value, m.params[i].error)

    reslist = [elem.value for elem in m.params]
    cov = np.array(m.covariance)
    jres.append(makeJfitres(reslist, cov, fitobs))


# Now saving the results:
j_val = {k : [] for k in jres[0]}
j_err = {k : [] for k in jres[0]}
for ires in jres:
    for iobs, ival in ires.items():
        j_val[iobs].append(ival.n)
        j_err[iobs].append(ival.s)

import json
outpath = f"output/fitres_N1e6_B{ctd_bins}_{ctl_bins}_{chi_bins}.json"
J = {"val" : j_val, "err": j_err, "bins" : b[0].tolist()}
with open(outpath, "w") as outfile:
    json.dump(J, outfile, indent=2)