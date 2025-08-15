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
from slophep.Fitting.PDFCoefSimple import PDFAngularCoef

from slophep.Predictions.Observables import BToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV.BToDstFF import BLPR
from slophep.Predictions.Math.BToVMathTools import calc_norm_j
import json

def makeJfitres(res, cov, order):
    j = uncert.correlated_values(res, cov)
    I1c = (4 - 6*j[0] + j[2] + 2*j[1])/3
    J = {order[i] : j[i] for i in range(len(j))}
    J['1c'] = I1c
    return J

wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.0,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.5
}
obs = BToDstEllNuPrediction("mu", "mu", BLPR)
obs.set_wc(wcoeffs)
ctd_bins = 3
ctl_bins = 3
chi_bins = 3
# Gets a histogram of the PDF with specified binning
h, b, angint, j_bins = obs.PDF_hist(10, ctd_bins, ctl_bins, chi_bins)

# Scale to approx expected yields
N = 6e5
h_norm = h/np.sum(h)
h_scale = N*h_norm
# Expect non-poisson errors + unfolding will increase uncertainty - set to 2*sigma_Poisson
h_scale_err = 2.0*np.sqrt(h_scale)
# h_scale_err = 0.06*h_scale

jres = []
fitobs = ['1s', '2s', '2c', '6s', '6c', 3, 4, 5, 7, 8, 9]
for ibin in range(len(h)):
    j_i = calc_norm_j(j_bins[ibin])
    h_data = h_scale[ibin]
    h_data_err = h_scale_err[ibin]
    jorder = angint["order"]

    # For simplicity set starting values at prediction
    init_vals = [j_i[k] for k in fitobs]
    PDF = PDFAngularCoef()
    PDF.setData(h_data, h_data_err)
    PDF.angint = angint["asarray"]
    reslist, cov = PDF.fit(init_vals)
    jres.append(makeJfitres(reslist, cov, fitobs))

    for i in range(len(reslist)):
        print(ibin, fitobs[i], j_i[fitobs[i]], init_vals[i], reslist[i], np.sqrt(cov[i][i]))

# Now saving the results:
j_val = {k : [] for k in jres[0]}
j_err = {k : [] for k in jres[0]}
for ires in jres:
    for iobs, ival in ires.items():
        j_val[iobs].append(ival.n)
        j_err[iobs].append(ival.s)


outpath = f"output/fitresNPT05_N6e5_B{ctd_bins}_{ctl_bins}_{chi_bins}.json"
J = {"val" : j_val, "err": j_err, "bins" : b[0].tolist()}
with open(outpath, "w") as outfile:
    json.dump(J, outfile, indent=2)