import numpy as np
import matplotlib.pyplot as plt

from slophep.Predictions.Observables import BuToD1stEllNuPrediction
from slophep.Predictions.FormFactorsBToDstst import BuToD1stFF
from slophep.Experimental.evtgen_correction import correction_weight_2body, correction_weight_3body, create_interpolated_weightfunc
from slophep.utils import setPlotParams

obs = BuToD1stEllNuPrediction("mu", "mu", BuToD1stFF.ISGW2)
gamma = obs.Gamma()

mvals = np.linspace(2., 5., 50)

setPlotParams()
# A D1'->2-body example
fig_2body, ax_2body = plt.subplots(1, 1)
# The correction_weight function obtains a weight for a particular mass
# This needs to integrate over q2 every function call which can be unoptimal
weights_2body = np.array([correction_weight_2body(im, obs, normscale=1./gamma) for im in mvals])
ax_2body.plot(mvals, weights_2body, label="Full")
# Can alternatively create a weight function that is constructed through interpolation
# this takes a bit of time to create the function, but the subsequent evaluations are faster
correction_func_2body = create_interpolated_weightfunc(correction_weight_2body, obs, 2., 5., Npoints=1000)
ax_2body.plot(mvals, correction_func_2body(mvals), 'k--', label="Interpolated")
ax_2body.legend()
ax_2body.set(xlabel = r"$m(D^{**})$ [GeV/$c$]", ylabel = "Weight (up to norm.)")
fig_2body.savefig("output/correction_weight_BuToD1pMuNu_D1p2body.png", dpi=100, bbox_inches="tight")


# A D1'->3-body example (D1p->D0PiPi):
fig_3body, ax_3body = plt.subplots(1, 1)
# This requires a bit more information for the extra phase-space factors
weights_3body = np.array([correction_weight_3body(im, obs, L=0, width=314.0e-3, daughters=["D0", "pi0", "pi0"], normscale=1./gamma) for im in mvals])
ax_3body.plot(mvals, weights_3body, label="Full")
# Can alternatively create a weight function that is constructed through interpolation
# this takes a bit of time to create the function, but the subsequent evaluations are faster
correction_func_3body = create_interpolated_weightfunc(
    correction_weight_3body, obs, 2., 5., Npoints=1000,
    weight_func_kwargs={"L" : 0, "width" : 314.0e-3, "daughters" : ["D0", "pi0", "pi0"]}
)
ax_3body.plot(mvals, correction_func_3body(mvals), 'k--', label="Interpolated")
ax_3body.legend()
ax_3body.set(xlabel = r"$m(D^{**})$ [GeV/$c$]", ylabel = "Weight (up to norm.)")
fig_3body.savefig("output/correction_weight_BuToD1pMuNu_D1p3body.png", dpi=100, bbox_inches="tight")
