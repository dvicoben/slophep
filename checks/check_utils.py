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


class ComparisonPlot:
    def __init__(self, ylabel: str):
        self.fig, self.ax = plt.subplots(1, 1)
        self.ax.set(xlabel=r"$q^2$ [GeV$^2$]", ylabel=ylabel)

    def add_slop_prediction(self, x: list[float], y: list[float], label: str):
        self.ax.plot(x, y, '-', label=label)
    
    def add_comparison_prediction(self, x: list[float], y: list[float], label: str):
        self.ax.plot(x, y, 'k--', label=label)
    
    def annotate(self, text: str, x: float, y: float):
        self.ax.text(x, y, text, size=8,
                    ha='left', va='top', transform=self.ax.transAxes)
    
    def makeplot(self):
        self.ax.legend(bbox_to_anchor=(1.01, 0.99), loc="upper left")
        return self.fig, self.ax
    
    def savefig(self, outpath: str, dpi: int = 150):
        self.fig.savefig(outpath, dpi=dpi, bbox_inches="tight")