import matplotlib.pyplot as plt
import numpy as np
import json
from bd2dstlnu.utils import setPlotParams

def readres(inpath: str, label: str = None):
    with open(inpath) as f:
        d = json.load(f)
    if label:
        d["label"] = label
    else:
        d["label"] = inpath
    return d


def makePlotHorizontal(obsdicts_list, labels, obslist, errorsonly, savepath):
    x_base = np.arange(1, len(obslist)+1, 1)
    plt.clf()
    plt.figure(figsize = (12, 4.5))
    maxabs = 0
    for i in range(len(obsdicts_list)):
        obs = obsdicts_list[i]["val"] if not errorsonly else np.zeros_like(obsdicts_list[i]["val"])
        err = obsdicts_list[i]["err"]
        x = x_base + i*0.5/(len(obsdicts_list)-1)
        if max(np.abs(obs)+np.abs(err)) > maxabs:
            maxabs = max(np.abs(obs)+np.abs(err))
        # print(obs)
        # plt.plot(x, obs, 'x', label=labels[i])
        plt.errorbar(x, obs, yerr=err, markersize=2, capsize=2, label=labels[i], fmt='o',)


    # plt.axhline(0.0, 0, 15, color='grey', alpha = 0.5, zorder=-5)
    tick_p = x_base + 0.25
    plt.gca().set_xticks(tick_p, obslist)
    
    ylim = [-maxabs-0.05, maxabs+0.05]
    plt.xlim(0.5, x_base[-1]+1)
    plt.ylim(ylim[0], ylim[1])

    for elem in tick_p:
        # plt.axhline(elem, color='grey', alpha = 0.5, linestyle='dashed')
        # plt.fill_between(ylim, elem-0.45, elem+0.45, zorder=-5, color='grey', alpha=0.1)
        plt.gca().axvspan(elem-0.45, elem+0.45, zorder=-5, color='grey', alpha=0.1)
    
    plt.legend(bbox_to_anchor=(1.01, 1.0))
    plt.tight_layout()
    plt.savefig(savepath, bbox_inches = 'tight', dpi=150)





setPlotParams()
results2221 = [
    "output/fitres_FFBGL2221_N2e6_B3_3_5.json",
    "output/fitres_FFBGL2221_N2e6_B4_4_5.json"
]
labels = [
    r"$N=2\times 10^6$, (3, 3, 5)",
    r"$N=2\times 10^6$, (4, 4, 5)",
]
res = [readres(i) for i in results2221]
makePlotHorizontal(res, labels, [r"$a_0$", r"$a_1$", r"$a_2$", r"$b_1$", r"$b_2$", r"$c_1$", r"$c_2$", r"$d_0$", r"$d_1$"], errorsonly=False, savepath="output/BGLfitres2221_compare.png")
makePlotHorizontal(res, labels, [r"$a_0$", r"$a_1$", r"$a_2$", r"$b_1$", r"$b_2$", r"$c_1$", r"$c_2$", r"$d_0$", r"$d_1$"], errorsonly=True, savepath="output/BGLfitres2221_compare_err.png")


results2221 = [
    "output/fitres_FFBGL1110_N2e6_B3_3_5.json",
    "output/fitres_FFBGL1110_N2e6_B4_4_5.json"
]
res = [readres(i) for i in results2221]
makePlotHorizontal(res, labels, [r"$a_0$", r"$a_1$", r"$b_1$", r"$c_1$", r"$d_0$"], errorsonly=False, savepath="output/BGLfitres1110_compare.png")
makePlotHorizontal(res, labels, [r"$a_0$", r"$a_1$", r"$b_1$", r"$c_1$", r"$d_0$"], errorsonly=True, savepath="output/BGLfitres1110_compare_err.png")