from slophep.Observables import BdToDstEllNuPrediction
from slophep.FormFactors import BdToDstFF

q2 = 5.0
obs = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.CLN())
# We save for later use here a dictionary of the default values
ffdefaults = obs.FF.define_userparams()


print("Observables at default values:")
print(obs.J(q2))
# Now change FF parameters to ones in HQET2 decfiles
hqet2 = {
    "RhoSq" : 1.122,
    "h_A1"  : 0.908,
    "R1"    : 1.270,
    "R2"    : 0.852,
    "R0"    : 1.15
}
obs.set_ff(hqet2)
print("Observables at HQET2 values:")
print(obs.J(q2))

# Now consider changing Coefficients to NP
# We are using CLN here which makes assumptions to get the tensor FFs, but this is simply illustrative
wcoeffs = {
    'CVL_bcmunumu': 0.0, 
    'CVR_bcmunumu': 0.5,
    'CSL_bcmunumu': 0.0,
    'CSR_bcmunumu': 0.0,
    'CT_bcmunumu': 0.0
}
obs.set_wc(wcoeffs)

obs.set_ff(ffdefaults)
print("Observables at default values and CVR=0.5:")
print(obs.J(q2))

obs.set_ff(hqet2)
print("Observables at HQET2 values and CVR=0.5:")
print(obs.J(q2))
