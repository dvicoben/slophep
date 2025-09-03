import numpy as np
import matplotlib.pyplot as plt
from slophep.Predictions.FormFactorBase import FormFactor

# Getting the SLOP predictions for BLPR:
def get_spectrum_slop(ff: FormFactor, 
                      qsq: list[float], 
                      method: str):
    res = {}
    for iq2 in qsq:
        ffmethod = getattr(ff, method)
        iff = ffmethod(iq2)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return {k : np.array(res[k]) for k in res}


def make_comparison_plot(sloppred: dict, 
                         otherpred: dict, 
                         qsq: list[float], 
                         otherpred_label: str,
                         ipreds: list[str],
                         ipredlabels: list[str] = [],
                         prefix: str = "",
                         savestr_fmt: str = "checks/check.png"):
    if len(ipredlabels) != len(ipreds):
        ipredlabels = ipreds
    for ipred, ilabel in zip(ipreds, ipredlabels):
        fig, ax = plt.subplots(1, 1)
        ax.plot(qsq, sloppred[ipred], 'b-', label="SLOP")
        ax.plot(qsq, otherpred[ipred], 'r--', label=otherpred_label)
        ax.set(xlabel=r"$q^2$", ylabel=ilabel, title=prefix)
        ax.legend()
        fig.savefig(savestr_fmt.format(prefix, ipred), 
                    bbox_inches = 'tight',
                    dpi=100)