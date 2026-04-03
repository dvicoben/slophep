import numpy as np
import matplotlib.pyplot as plt

from slophep.Predictions.Observables import BuToD1stEllNuPrediction
from slophep.Predictions.FormFactorsBToDstst import BuToD1stFF
from slophep.Experimental.evtgen_correction import correction_weight, create_interpolated_weightfunc
from slophep.utils import setPlotParams

obs = BuToD1stEllNuPrediction("mu", "mu", BuToD1stFF.ISGW2)
gamma = obs.Gamma()

mvals = np.linspace(2., 5., 50)
setPlotParams()
fig, ax = plt.subplots(1, 1)
# The correction_weight function obtains a weight for a particular mass
# This needs to integrate over q2 every function call which can be unoptimal
weights = np.array([correction_weight(im, obs, 1./gamma) for im in mvals])
ax.plot(mvals, weights, label="Full")
# Can alternatively create a weight function that is constructed through interpolation
# this takes a bit of time to create the function, but the subsequent evaluations are faster
correction_func = create_interpolated_weightfunc(obs, 2., 5., Npoints=1000)
ax.plot(mvals, correction_func(mvals), 'k--', label="Interpolated")

ax.set(xlabel = r"$m(D^{**})$ [GeV/$c$]", ylabel = "Weight (up to norm.)")
fig.savefig("output/correction_weight_BuToD1pMuNu.png", dpi=100, bbox_inches="tight")