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
        

def get_additional_spectrum_BGL(qsq, ff):
    res = {"z" : []}
    for iq2 in qsq:
        w = ff.w(iq2)
        z = (np.sqrt(w+1) - np.sqrt(2))/(np.sqrt(w+1) + np.sqrt(2))
        res["z"].append(z)
        Mb = ff.internalparams["Mb"]
        Mc = ff.internalparams["Mc"]
        for ielem in ["f", "g", "P1"]:
            if f"Bl{ielem}" not in res:
                res[f"Bl{ielem}"] = []
            res[f"Bl{ielem}"].append(ff.blaschke(ff.internalparams[f"BcStates{ielem}"], z, Mb, Mc))

        outer_fcn = ff.outer_fcns(z)
        for elem in outer_fcn:
            if elem not in res:
                res[elem] = []
            res[elem].append(outer_fcn[elem])

    res = {k : np.array(res[k]) for k in res}
    return res