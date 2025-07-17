from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.Fluctuate import FluctuatePrediction
from bd2dstlnu.Predictions.BToDstFFHPQCD import HPQCD
from bd2dstlnu.Predictions.BToDstFFBSZ import BSZ
import json
import numpy as np
import matplotlib.pyplot as plt

wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.0,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.0
}

obs_hpqcd = BToDstEllNuPrediction("mu", "mu", HPQCD)
obs_hpqcd.set_wc(wcoeffs)
obs_bsz = BToDstEllNuPrediction("mu", "mu", BSZ)
obs_bsz.set_wc(wcoeffs)

f_hpqcd = FluctuatePrediction(obs_hpqcd)
fconfig_hpqcd = "data/FF_HPQCD_COV_arXiv230403137.json"
f_hpqcd.set_params_from_configfile(fconfig_hpqcd)
f_hpqcd.fluctuate(1000)
jerr_hpqcd = f_hpqcd.get_error("J", [5.0])


f_bsz = FluctuatePrediction(obs_bsz)
fconfig_bsz = "data/FF_BSZ_COV_arXiv181100983.json"
f_bsz.set_params_from_configfile(fconfig_bsz)
f_bsz.fluctuate(1000)
jerr_bsz = f_bsz.get_error("J", [5.0])

print(obs_hpqcd.J(5.0))
print(jerr_hpqcd)
print(obs_bsz.J(5.0))
print(jerr_bsz)


# Now th FFs themselves
ff_hpqcd = FluctuatePrediction(HPQCD())
ff_hpqcd.set_params_from_configfile(fconfig_hpqcd)
ff_hpqcd.fluctuate(1000)
fferr_hpqcd = ff_hpqcd.get_error("get_ff", [5.0])


ff_bsz = FluctuatePrediction(BSZ())
ff_bsz.set_params_from_configfile(fconfig_bsz)
ff_bsz.fluctuate(1000)
fferr_bsz = ff_bsz.get_error("get_ff", [5.0])


hpqcd = HPQCD()
bsz = BSZ()

def get_spectrum(ff, fluct, npoints = 50):
    qsq = np.linspace(obs_hpqcd.q2min, obs_hpqcd.q2max, npoints)
    res = {}
    for iq2 in qsq:
        qcd_f = ff.get_ff(iq2)
        qcd_f_err = fluct.get_error("get_ff", [iq2])
        for iff in qcd_f:
            if iff not in res:
                res[iff] = {"val" : [], "lo" : [], "hi" : []}
            res[iff]["val"].append(qcd_f[iff])
            res[iff]["lo"].append(qcd_f_err[iff][0])
            res[iff]["hi"].append(qcd_f_err[iff][1])
    return qsq, res


qsq, h_res = get_spectrum(hpqcd, ff_hpqcd)
qsq, b_res = get_spectrum(bsz, ff_bsz)
for iff, ires in h_res.items():
    plt.plot(qsq, ires["val"], 'b-')
    plt.fill_between(qsq, ires["lo"], ires["hi"], alpha=0.15, color='b', label="HPQCD")
    if iff in b_res:
        plt.plot(qsq, b_res[iff]["val"], 'r-')
        plt.fill_between(qsq, b_res[iff]["lo"], b_res[iff]["hi"], alpha=0.15, color='r', label="BSZ")
    plt.legend()
    plt.savefig(f"output/load_fit_comp_{iff}.png", bbox_inches='tight')
    plt.clf()
    plt.close()