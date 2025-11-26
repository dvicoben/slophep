"""
Simple script to test FF fits using BToDstObs prediction as fit model

NOTE: At the moment fits are quite slow as:
    - Method (as implemented) recomputes the angular integrals at every step 
      (can be reduced to compute only once, since these are always the same)
    - Need to recompute and integrate J_i at every fit step (cannot really be remedied)
"""

import numpy as np
from slophep.Fitting.PDFBase import PDFBase
from slophep.Fitting.Parameter import Parameter, ParameterManager
from slophep.Fitting.CostFuncs import Chi2Uncorrelated, Chi2Correlated
from slophep.Fitting.Fitter import imFitter
from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV.BdToDstFF import BGL, BGL_Hammer
import matplotlib.pyplot as plt
from slophep.utils import setPlotParams

obs = BdToDstEllNuPrediction("mu", "mu", BGL)
print("Generated Values:")
print(obs.FF.ffpar)
h, b = obs.dGdq2_hist(10)
bc = 0.5*(b[1:] + b[:-1])
bw = 0.5*(b[1:] - b[:-1])
# Gets a histogram of the PDF with specified binning
N = 1e6
h_norm = h/np.sum(h)
h_data = N*h_norm
# Expect non-poisson errors + unfolding will increase uncertainty - set to 2*sigma_Poisson
h_data_err = 2.0*np.sqrt(h_data)



class BGLPDF(PDFBase):
    def __init__(self, name: str, predobj: BdToDstEllNuPrediction):
        super().__init__(name, predobj)
        # self._params_list = ["a0", "a1", "b1", "c1", "d0"]
        self._params_list = ["a0", "a1", "b1", "c1", "d0"]
        self.bc = None
        self.b = None

    def pdf(self):
        self.predobj.set_ff({k : self.getVal(k) for k in self.params_list})
        # x = [self.predobj.dGdq2(iq2) for iq2 in self.bc]
        x, _ = self.predobj.dGdq2_hist(self.b)
        return np.array(x)


pred = BdToDstEllNuPrediction("mu", "mu", BGL)
pdf = BGLPDF("BGLPDF", pred)
pdf.bc = 0.5*(b[1:] + b[:-1])
pdf.b = b

param_manager = ParameterManager()
for ielem in pdf.params_list:
    ipar = Parameter(ielem, obs.FF.ffpar[ielem])
    if ielem == "b0":
        ipar.fix()
    param_manager.addParam(ipar)

pdf.link_params({k : k for k in pdf.params_list})
chi2 = Chi2Uncorrelated(pdf, param_manager, h_data, h_data_err, True)
m = imFitter(chi2)
reslist, err = m.fit(False)
for i in range(len(pdf.params_list)):
    print(pdf.params_list[i], reslist[i], err[i], obs.FF.ffpar[pdf.params_list[i]])



setPlotParams()
fig, ax = plt.subplots()
h_fitres, b = pdf.predobj.dGdq2_hist(b)
h_fitres_scaled = N*h_fitres/np.sum(h_fitres)
ax.errorbar(bc, h_data, h_data_err, bw, "ko", capsize=2, markersize=2, label="Toy")
ax.stairs(h_fitres_scaled, b, color="b", label="Fit")
ax.set(xlim=(b[0], b[-1]), xlabel=r"$q^2$ [GeV$^2$]")
ax.legend()
fig.savefig("output/fitBGL_plot.png", bbox_inches="tight", dpi=100)