import numpy as np

from slophep.Generator import MCGenerator
from slophep.Predictions.Observables import BToPEllNuPrediction, BToVEllNuPrediction


class MCGeneratorBToV(MCGenerator):
    def __init__(self, obs: BToVEllNuPrediction, seed: int = None):
        super().__init__(obs, seed)
    
    def get_random_point(self) -> list[float]:        
        q2min, q2max = self.obs.q2min, self.obs.q2max
        p = [
            self.rng.uniform(q2min, q2max),    # q2
            self.rng.uniform(-1, 1),           # ctx
            self.rng.uniform(-1, 1),           # ctl
            self.rng.uniform(-np.pi, np.pi)    # chi
        ]
        return p
    

class MCGeneratorBToP(MCGenerator):
    def __init__(self, obs: BToPEllNuPrediction, seed: int = None):
        super().__init__(obs, seed)
    
    def get_random_point(self) -> list[float]:        
        q2min, q2max = self.obs.q2min, self.obs.q2max
        p = [
            self.rng.uniform(q2min, q2max),    # q2
            self.rng.uniform(-1, 1)            # ctl
        ]
        return p