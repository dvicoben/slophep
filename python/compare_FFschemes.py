from bd2dstlnu.Predictions.Observables.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.utils import setPlotParams
from bd2dstlnu.Predictions.FormFactorsBToDst import CLN, BGL, BLPR, BSZ, HPQCD

import matplotlib.pyplot as plt
        
setPlotParams()

wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.0,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.0
}

obs_cln = BToDstEllNuPrediction("mu", "mu", CLN)
obs_blpr = BToDstEllNuPrediction("mu", "mu", BLPR)
obs_bgl = BToDstEllNuPrediction("mu", "mu", BGL)
obs_bsz = BToDstEllNuPrediction("mu", "mu", BSZ)
obs_hpqcd = BToDstEllNuPrediction("mu", "mu", HPQCD)

obs_cln.set_wc(wcoeffs)
obs_blpr.set_wc(wcoeffs)
obs_bgl.set_wc(wcoeffs)
obs_bsz.set_wc(wcoeffs)
obs_hpqcd.set_wc(wcoeffs)

print("Test at q2=5.0:")
print(obs_cln.J(5.0))
print(obs_blpr.J(5.0))
print(obs_bgl.J(5.0))
print(obs_bsz.J(5.0))
print(obs_hpqcd.J(5.0))

obs = "FL"
p = obs_cln.plot_obs_prediction(obs, label="CLN")
p = obs_blpr.plot_obs_prediction(obs, label="BLPR", plot=p)
p = obs_bgl.plot_obs_prediction(obs, label="BGL", plot=p)
p = obs_bsz.plot_obs_prediction(obs, label="BSZ", plot=p)
p = obs_hpqcd.plot_obs_prediction(obs, label="HPQCD", plot=p)
p[1].set(ylabel=r"$F_L$")
p[1].legend()
p[0].savefig(f"output/{obs}_prediction.png", bbox_inches="tight")
p[0].show()