from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.SamplingFluctuate import SamplingHelper
from bd2dstlnu.Predictions.BToDstFFHPQCD import HPQCD
from bd2dstlnu.Predictions.BToDstFFBSZ import BSZ

from bd2dstlnu.utils import setPlotParams
import numpy as np
import matplotlib.pyplot as plt

Nfluct = 5000
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

f_hpqcd = SamplingHelper(obs_hpqcd)
fconfig_hpqcd = "data/FF_HPQCD_COV_arXiv230403137.json"
f_hpqcd.set_params_from_configfile(fconfig_hpqcd)
f_hpqcd.fluctuate(Nfluct)
jerr_hpqcd = f_hpqcd.get_error("J", [5.0])


f_bsz = SamplingHelper(obs_bsz)
fconfig_bsz = "data/FF_BSZ_COV_arXiv181100983.json"
f_bsz.set_params_from_configfile(fconfig_bsz)
f_bsz.fluctuate(Nfluct)
jerr_bsz = f_bsz.get_error("J", [5.0])

print(obs_hpqcd.J(5.0))
print(jerr_hpqcd)
print(obs_bsz.J(5.0))
print(jerr_bsz)


# Now th FFs themselves
ff_hpqcd = SamplingHelper(HPQCD())
ff_hpqcd.set_params_from_configfile(fconfig_hpqcd)
ff_hpqcd.fluctuate(Nfluct)
fferr_hpqcd = ff_hpqcd.get_error("get_ff", [5.0])


ff_bsz = SamplingHelper(BSZ())
ff_bsz.set_params_from_configfile(fconfig_bsz)
ff_bsz.fluctuate(Nfluct)
fferr_bsz = ff_bsz.get_error("get_ff", [5.0])


hpqcd = HPQCD()
bsz = BSZ()


def get_spectrum_dict(qsq, obs, attr, fluct):
    res = {}
    for iq2 in qsq:
        o = getattr(obs, attr)(iq2)
        o_err = fluct.get_error(attr, [iq2])
        for iobs in o:
            if iobs not in res:
                res[iobs] = {"val" : [], "lo" : [], "hi" : []}
            res[iobs]["val"].append(o[iobs])
            res[iobs]["lo"].append(o_err[iobs][0])
            res[iobs]["hi"].append(o_err[iobs][1])
    return res

def plot_spectrum_dict(qsq, obslist, res_list, label_list):
    
    for iobs in obslist:
        for res, label in zip(res_list, label_list):
            if iobs not in res:
                continue
            ires = res[iobs]
            lines = plt.plot(qsq, ires["val"])
            fill = plt.fill_between(qsq, ires["lo"], ires["hi"], color=lines[0].get_color(), alpha=0.15, label=label)
        
        plt.xlabel(r"$q^2$")
        plt.ylabel(iobs)
        plt.legend()
        plt.savefig(f"output/plot_q2sepctrum_{iobs}.png", bbox_inches='tight')
        plt.clf()
        plt.close()



npoints = 100
qsq = np.linspace(obs_hpqcd.q2min, obs_hpqcd.q2max, npoints)

setPlotParams()

hpqcd_ff_res = get_spectrum_dict(qsq, hpqcd, "get_ff", ff_hpqcd)
bsz_ff_res = get_spectrum_dict(qsq, bsz, "get_ff", ff_bsz)
plot_spectrum_dict(qsq, [ielem for ielem in hpqcd_ff_res], 
                   [hpqcd_ff_res, bsz_ff_res], ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"])

hpqcd_J_res = get_spectrum_dict(qsq, obs_hpqcd, "J", f_hpqcd)
bsz_J_res = get_spectrum_dict(qsq, obs_bsz, "J", f_bsz)
plot_spectrum_dict(qsq, [ielem for ielem in hpqcd_J_res], 
                   [hpqcd_J_res, bsz_J_res], ["HPQCD arXiv:2304.03137", "BSZ arXiv:1811.00983"])