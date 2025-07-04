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

from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.BToDstFFBLPR import BLPR
from bd2dstlnu.Predictions.BToDstMathTools import calc_norm_j


def makeJfitres(res, cov, order):
    j = uncert.correlated_values(res, cov)
    I1c = (4 - 6*j[0] + j[2] + 2*j[1])/3
    J = {order[i] : j[i] for i in range(len(j))}
    J['1c'] = I1c
    return J

# wcoeffs = {
#     'CVL_bcmunumu': 0.0, 
#     'CVR_bcmunumu': 0.0,
#     'CSL_bcmunumu': 0.0,
#     'CSR_bcmunumu': 0.0,
#     'CT_bcmunumu': 0.0
# }

obs = BToDstEllNuPrediction("mu", "mu", BLPR)
# Gets a histogram of the PDF with specified binning
h, b, angint, j_bins = obs.PDF_hist(10, 4, 4, 5)

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


# Now plotting the results:

def setPlotParams(params = {}):
    """Set of parameters for plots in matplotlib
    """
    if len(params) > 0:
        plt.rcParams.update(params)
    else:
        font = {'family' : 'cmr10', #Or Times New Roman
            'size'   : 14}
        plt.rc('font', **font)
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams["axes.formatter.use_mathtext"] = True
        plt.rcParams.update({'axes.labelsize': 14,
                            'legend.fontsize': 10,
                            'xtick.labelsize': 12,
                            'ytick.labelsize': 12,
                            'figure.figsize': [8, 8/1.618]})

setPlotParams()


def plotresults(j_dicts: list, q2_bins: list):
    j_val = {k : [] for k in j_dicts[0]}
    j_err = {k : [] for k in j_dicts[0]}
    for ires in j_dicts:
        for iobs, ival in ires.items():
            j_val[iobs].append(ival.n)
            j_err[iobs].append(ival.s)

    bc = 0.5*(q2_bins[1:] + q2_bins[:-1])
    bw = 0.5*(q2_bins[1:] - q2_bins[:-1])

    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axs = axes.flat
    ax_idx = 0
    kwargs = {"fmt" : 'o', "capsize" : 2, "markersize" : 2, "color" : 'k'}

    for iobs in j_val:
        ax = axs[ax_idx]
        ax.errorbar(bc, j_val[iobs], yerr=j_err[iobs], xerr=bw, **kwargs)
        ax.set(xlabel = r"$q^2$ [GeV$^2$]", ylabel=r"$J_{" + str(iobs) + r"}$",
               xlim=(q2_bins[0], q2_bins[-1]))
        ax_idx += 1

    fig.tight_layout()
    return fig, axes


p = plotresults(jres, b[0])
p[0].savefig("output/Jfitres.png", bbox_inches='tight')