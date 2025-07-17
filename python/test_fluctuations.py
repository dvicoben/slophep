from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.Fluctuate import FluctuatePrediction
from bd2dstlnu.Predictions.BToDstFFHPQCD import HPQCD
from bd2dstlnu.Predictions.BToDstFFBSZ import BSZ
import json

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
f_hpqcd.fluctuate(500)
jerr_hpqcd = f_hpqcd.get_error("J", [5.0])


f_bsz = FluctuatePrediction(obs_bsz)
fconfig_bsz = "data/FF_BSZ_COV_arXiv181100983.json"
f_bsz.set_params_from_configfile(fconfig_bsz)
f_bsz.fluctuate(500)
jerr_bsz = f_bsz.get_error("J", [5.0])

print(obs_hpqcd.J(5.0))
print(jerr_hpqcd)
print(obs_bsz.J(5.0))
print(jerr_bsz)


# Now th FFs themselves
ff_hpqcd = FluctuatePrediction(HPQCD())
ff_hpqcd.set_params_from_configfile(fconfig_hpqcd)
ff_hpqcd.fluctuate(500)
fferr_hpqcd = ff_hpqcd.get_error("get_ff", [5.0])


ff_bsz = FluctuatePrediction(BSZ())
ff_bsz.set_params_from_configfile(fconfig_bsz)
ff_bsz.fluctuate(500)
fferr_bsz = ff_bsz.get_error("get_ff", [5.0])