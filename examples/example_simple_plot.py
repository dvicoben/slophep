from slophep.Observables import BdToDstEllNuPrediction
from slophep.FormFactors import BdToDstFF


obs_cln = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.CLN())
# We save for later use here a dictionary of the default values
ffdefaults = obs_cln.FF.define_userparams()

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
p[0].savefig("output/example_simple_plot.png", bbox_inches = "tight")
p[0].show()