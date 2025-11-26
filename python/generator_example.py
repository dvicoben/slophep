from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV import BdToDstFF
from slophep.Generator import MCGeneratorBToV
from slophep.Generator.root_interface import generate_to_root

obs = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.BLPR)
gen = MCGeneratorBToV(obs, 50)
generate_to_root(gen, 4000, "ignore/MC_BLPR_test.root", printprog=10)