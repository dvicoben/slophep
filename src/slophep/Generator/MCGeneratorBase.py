from slophep.Predictions.Observables import ObservableBase
import random

# from typing import Any
# PredictionObj = Any # in principle BToVEllNuPrediction or BToPEllNuPrediction

class MCGenerator:
    def __init__(self, obs: ObservableBase, seed: int = None):
        self._obs = obs
        self._maxBF = -1
        self._maxprob = 0
        # self._point_max_attempts = None
        self._seed = seed
        self._rng = random.Random(seed)
        self._coords = []

    @property
    def obs(self) -> ObservableBase:
        """The HbToHEllNuPrediction object"""
        return self._obs
    @property
    def maxBF(self) -> float:
        """Maximum BF to use in accept/reject"""
        return self._maxBF
    @property
    def maxprob(self) -> float:
        """Maximum probability to use in accept/reject"""
        return self._maxprob
    @property
    def seed(self) -> int:
        """RNG seed"""
        return self._seed
    @property
    def rng(self) -> random.Random:
        """The random number generator"""
        return self._rng
    @property
    def coords(self) -> list[str]:
        """Names of each coordinate to fully specify a point"""
        return self._coords
    
    def get_max_BF(self):
        """Finds maximum BF in q2 for usage in rejection method"""
        q2min, q2max = self.obs.q2min, self.obs.q2max
        maxBF = 0
        for iq2 in range(1000):
            q2 = q2min + iq2*(q2max - q2min)/float(1000-1)
            bf = self.obs.dBRdq2(q2)
            maxBF = max(bf, maxBF)
        self._maxBF = maxBF
        return maxBF
    
    def get_random_point(self) -> list[float]:
        """Get a randomly generated point (uniform). Implemented in derived class as different modes have different angular variables.

        Returns
        -------
        list[float]
            Point, coordinates in order of self.coords

        Raises
        ------
        NotImplementedError
            Must be implemented in dervied class
        """
        raise NotImplementedError("get_random_point method implemented in derived class!")

    def generate_point(self) -> list[float]:
        """Generate a point according to Prediction object PDF using accept/reject method

        Returns
        -------
        list[float]
            Point, coordinates in order of self.coords

        Raises
        ------
        ValueError
            If maxBF has not been run, cannot use rejection method, must run maxBF before generating points.
        """
        if self.maxBF < 0:
            raise ValueError(f"Max BF for rejection method is {self.maxBF}, unphysical. Remember to run get_max_BF.")
        point = []
        while len(point) < 1:
            ptrial = self.get_random_point()
            pbf = self.obs.dBRdq2(ptrial[0])
            pq2 = self.rng.uniform(0, 1)
            if pq2 > pbf:
                continue
            p = self.obs.PDF_norm(*ptrial)
            if p > self.maxprob:
                self._maxprob = p
            if p <=0: 
                continue
            prob = self.rng.uniform(0, 1.66)
            if prob < p:
                point = ptrial
        return point

    def generateN(self, N: int) -> list:
        """Generate N points according to prediction object PDF

        Parameters
        ----------
        N : int
            Number if points to generate

        Returns
        -------
        list
            List of points
        """
        npoints = []
        for i in range(N):
            npoints.append(self.generate_point())
        return npoints