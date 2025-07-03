from bd2dstlnu.Predictions.BToDstObs import BToDstEllNuPrediction
from bd2dstlnu.Predictions.BToDstFFBLPR import BLPR

# wcoeffs = {
#     'CVL_bcmunumu': 0.0, 
#     'CVR_bcmunumu': 0.0,
#     'CSL_bcmunumu': 0.0,
#     'CSR_bcmunumu': 0.0,
#     'CT_bcmunumu': 0.0
# }

obs = BToDstEllNuPrediction("mu", "mu", BLPR)
# Gets a histogram of the PDF with specified binning
h, b, angint = obs.PDF_hist(10, 4, 4, 5)