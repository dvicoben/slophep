import numpy as np
import matplotlib.pyplot as plt
import pypmc

from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV import BdToDstFF

from slophep.Predictions.MCGenerator import MCGenerator

pred = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.CLN)
norm = pred.dGdq2_bin(pred.q2min, pred.q2max)


lower_bnds = [pred.q2min, -1, -1, -np.pi]
upper_bnds = [pred.q2max, 1, 1, np.pi]

gen = MCGenerator(230)
gen.setPDF(pred.PDF, [(ilo, ihi) for ilo, ihi in zip(lower_bnds, upper_bnds)])

N       = 60000  # Number of samples that shall be returned
stride  = 5      # Stride, i.e., the number by which the actual amount of samples shall be thinned to return N samples
preruns = 8      # Number of preruns
pre_N   = 3000   # Number of samples per prerun

gen.init_sampler()
gen.prerun(preruns, pre_N)
parameter_samples = gen.sample(N, stride)


fig, axs = plt.subplots(1, 4, figsize=(15, 5))
Nbins = 30
ax = axs.flat
predfuncs = [
    lambda x : pred.dGdq2(x)/norm,
    pred.dGdctv_norm,
    pred.dGdctl_norm,
    pred.dGdchi_norm,
]

for i in range(len(axs.flat)):
    iax = axs[i]
    varidx = i
    dat = parameter_samples.T[varidx]
    ibnds = (lower_bnds[varidx], upper_bnds[varidx])
    h, b = np.histogram(dat, Nbins, range=ibnds, density=True)
    hraw, _ = np.histogram(dat, Nbins, range=ibnds)
    hrawerr = np.sqrt(hraw)
    hrawerr *= h/hraw
    bc = 0.5*(b[1:] + b[:-1])
    bw = 0.5*(b[1:] - b[:-1])
    ax[i].errorbar(bc, h, yerr=hrawerr, xerr = bw, fmt="ko", markersize=3, capsize=2)
    predpoints = np.linspace(*ibnds, 50)
    predv = [predfuncs[i](ip) for ip in predpoints]
    ax[i].plot(predpoints, predv, "r-")
    ax[i].set_ylim(bottom=0)
    ax[i].set(xlim=ibnds)
    # ax[i].hist(parameter_samples.T[varidx], bins=Nbins, range=(lower_bnds[varidx], upper_bnds[varidx]))
