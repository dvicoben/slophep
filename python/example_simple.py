from slophep.Predictions.Observables import BToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV import BToDstFF


obs_cln = BToDstEllNuPrediction("mu", "mu", BToDstFF.CLN2)
print(obs_cln.FF.ffpar) # print the default FF parameters
ffdefaults = obs_cln.FF.ffpar.copy()

# Plot out the FL prediction
p = obs_cln.plot_obs_prediction("FL", label="SM, Hammer defaults")
p[1].set(ylabel = r"$F_L$")

# Now change FF parameters to ones in HQET2 decfiles
hqet2 = {
    "RhoSq" : 1.122,
    "h_A1"  : 0.908,
    "R1"    : 1.270,
    "R2"    : 0.852,
    "R0"    : 1.15
}
obs_cln.set_ff(hqet2)
print(obs_cln.FF.ffpar) # print FF parameters again - should have changed!
p = obs_cln.plot_obs_prediction("FL", label="SM, HQET2", plot=p)

# Now consider changing Coefficients to NP
# We are using CLN here which makes assumptions to get the tensor FFs, but this is simply illustrative
wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.5,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.0
}
obs_cln.set_wc(wcoeffs)

# Lets plot now NP for both sets of ff parameters
obs_cln.set_ff(ffdefaults)
p = obs_cln.plot_obs_prediction("FL", label="NP, Hammer defaults", plot=p)
obs_cln.set_ff(hqet2)
p = obs_cln.plot_obs_prediction("FL", label="NP, HQET2", plot=p)
p[1].legend()
p[0].show()