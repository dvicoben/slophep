import matplotlib.pyplot as plt
import numpy as np
import json
from slophep.utils import setPlotParams

def readres(inpath: str, label: str = None):
    with open(inpath) as f:
        d = json.load(f)
    if label:
        d["label"] = label
    else:
        d["label"] = inpath
    return d


def plotresults(j_val: dict, j_err: dict, q2_bins: list):
    q2_bins = np.array(q2_bins)
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


def plotresult_multiple(jres: list):
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axs = axes.flat
    kwargs = {"fmt" : 'o', "capsize" : 2, "markersize" : 2}

    for elem in jres:
        j_val = elem["val"]
        j_err = elem["err"]
        q2_bins = np.array(elem["bins"])
        bc = 0.5*(q2_bins[1:] + q2_bins[:-1])
        bw = 0.5*(q2_bins[1:] - q2_bins[:-1])
        
        ax_idx = 0
        for iobs in j_val:
            ax = axs[ax_idx]
            ax.errorbar(bc, j_val[iobs], yerr=j_err[iobs], xerr=bw, label=elem["label"], **kwargs)
            ax_idx += 1
            ax.set(xlabel = r"$q^2$ [GeV$^2$]", ylabel=r"$J_{" + str(iobs) + r"}$")
    
    axs[0].legend()
    fig.tight_layout()

    return fig, axes


def plotresult_multiple_errors_only(jres: list):
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axs = axes.flat
    kwargs = {"fmt" : 'o', "capsize" : 3, "markersize" : 2}

    for elem in jres:
        j_val = elem["val"]
        j_err = elem["err"]
        q2_bins = np.array(elem["bins"])
        bc = 0.5*(q2_bins[1:] + q2_bins[:-1])
        bw = 0.5*(q2_bins[1:] - q2_bins[:-1])
        
        ax_idx = 0
        for iobs in j_val:
            ax = axs[ax_idx]
            ax.errorbar(bc, np.zeros_like(j_val[iobs]), yerr=j_err[iobs], xerr=bw, label=elem["label"], **kwargs)
            ax_idx += 1
            ax.set(xlabel = r"$q^2$ [GeV$^2$]", ylabel=r"$J_{" + str(iobs) + r"}$")
    
    axs[0].legend()
    fig.tight_layout()

    return fig, axes


setPlotParams()

# r = readres("output/fitres_N1e6_B3_3_5.json")
# p = plotresults(r["val"], r["err"], r["bins"])
# p[0].savefig("output/Jfitres_0.png", bbox_inches='tight')

# results = [
#     "output/fitres_N3e5_B3_3_5.json",
#     "output/fitres_N3e5_B4_4_5.json",
#     "output/fitres_N3e5_B5_5_5.json",
#     "output/fitres_N1e6_B3_3_5.json",
#     "output/fitres_N1e6_B4_4_5.json",
#     "output/fitres_N1e6_B5_5_5.json"
# ]
# labels = [
#     r"$N=3\times 10^5$, (3, 3, 5)",
#     r"$N=3\times 10^5$, (4, 4, 5)",
#     r"$N=3\times 10^5$, (5, 5, 5)",
#     r"$N=1\times 10^6$, (3, 3, 5)",
#     r"$N=1\times 10^6$, (4, 4, 5)",
#     r"$N=1\times 10^6$, (5, 5, 5)"
# ]
results = [
    "output/fitres_N6e5_B3_3_3.json",
    "output/fitres_N6e5_B3_3_5.json",
    "output/fitres_N6e5_B4_4_5.json",
    "output/fitres_N2e6_B3_3_3.json",
    "output/fitres_N2e6_B3_3_5.json",
    "output/fitres_N2e6_B4_4_5.json",
]
labels = [
    r"$N=6\times 10^5$, (3, 3, 3)",
    r"$N=6\times 10^5$, (3, 3, 5)",
    r"$N=6\times 10^5$, (4, 4, 5)",
    r"$N=2\times 10^6$, (3, 3, 3)",
    r"$N=2\times 10^6$, (3, 3, 5)",
    r"$N=2\times 10^6$, (4, 4, 5)",
]
jres = [readres(ip, ilabel) for ip, ilabel in zip(results, labels)]
p = plotresult_multiple(jres)
p[0].savefig("output/Jfitres_compare.png", bbox_inches='tight')

perr = plotresult_multiple_errors_only(jres)
perr[0].savefig("output/Jfitres_compare_err.png", bbox_inches='tight')


results = [
    "output/fitresNPVRj05_N6e5_B3_3_3.json",
    "output/fitresNPVRj05_N6e5_B3_3_5.json",
    "output/fitresNPVRj05_N6e5_B4_4_5.json",
    "output/fitresNPVRj05_N2e6_B3_3_3.json",
    "output/fitresNPVRj05_N2e6_B3_3_5.json",
    "output/fitresNPVRj05_N2e6_B4_4_5.json",
]
jres = [readres(ip, ilabel) for ip, ilabel in zip(results, labels)]
p = plotresult_multiple(jres)
p[0].savefig("output/Jfitres_NPVR_compare.png", bbox_inches='tight')

perr = plotresult_multiple_errors_only(jres)
perr[0].savefig("output/Jfitres_NPVR_compare_err.png", bbox_inches='tight')


results = [
    "output/fitresNPT05_N6e5_B3_3_3.json",
    "output/fitresNPT05_N6e5_B3_3_5.json",
    "output/fitresNPT05_N6e5_B4_4_5.json",
    "output/fitresNPT05_N2e6_B3_3_3.json",
    "output/fitresNPT05_N2e6_B3_3_5.json",
    "output/fitresNPT05_N2e6_B4_4_5.json",
]
jres = [readres(ip, ilabel) for ip, ilabel in zip(results, labels)]
p = plotresult_multiple(jres)
p[0].savefig("output/Jfitres_NPT_compare.png", bbox_inches='tight')

perr = plotresult_multiple_errors_only(jres)
perr[0].savefig("output/Jfitres_NPT_compare_err.png", bbox_inches='tight')