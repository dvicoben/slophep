from typing import Callable
import numpy as np
import pypmc

import logging
logger = logging.getLogger(__name__)

SamplingFunction = Callable[[*tuple[float, ...]], float]

class MCGenerator:
    def __init__(self, seed: int = None):
        self._seed    : int                          = seed
        self._pdf     : SamplingFunction             = None
        self._bounds  : list[tuple[float, float]]    = []
        self._rng     : np.random.RandomState        = np.random.mtrand if seed is None else np.random.mtrand.RandomState(seed)
        self._sampler : pypmc.sampler.markov_chain.AdaptiveMarkovChain = None

    @property
    def pdf(self) -> SamplingFunction:
        return self._pdf
    @property
    def bounds(self) -> list[tuple[float, float]]:
        return self._bounds
    @property
    def seed(self) -> int:
        return self._seed
    @property
    def rng(self) -> np.random.RandomState:
        return self._rng
    @property
    def sampler(self) -> pypmc.sampler.markov_chain.AdaptiveMarkovChain:
        return self._sampler

    def logpdf(self, x: list[float]) -> float:
        return np.log(self.pdf(*x))
    
    def setPDF(self, pdf: SamplingFunction, bounds: list[tuple[float, float]]) -> None:
        """Set the pdf to sample

        Parameters
        ----------
        pdf : SamplingFunction
            Callable to sample, should have float arguments pdf(x: float, y: float, z: float)
        bounds : list[tuple[float, float]]
            Bounds for arguments of pdf, in positional order
        """
        self._pdf    = pdf
        self._bounds = bounds

    def init_sampler(self, start_point: list[float] = None, cov_scale: float = 0.1):
        bnd_lo = np.array([ibnd[0] for ibnd in self.bounds])
        bnd_up = np.array([ibnd[1] for ibnd in self.bounds])
        bnd = pypmc.tools.indicator.hyperrectangle(bnd_lo, bnd_up)
        log_target = pypmc.tools.indicator.merge_function_with_indicator(self.logpdf, bnd, -np.inf)
        sigma = np.diag([np.square(bnd_lo[k] - bnd_up[k]) / 12 * cov_scale for k in range(len(bnd_lo))])
        log_proposal = pypmc.density.gauss.LocalGauss(sigma)

        if start_point is None:
            u = np.array([self.rng.uniform(0.0, 1.0) for j in range(0, len(bnd_lo))])
            ubar = 1.0 - u
            start_point = ubar*bnd_up + u*bnd_lo
        else:
            start_point = np.array(start_point)
        
        self._sampler = pypmc.sampler.markov_chain.AdaptiveMarkovChain(log_target, log_proposal, start_point, save_target_values=True, rng=self.rng)

    def prerun(self, preruns: int, preN: int) -> None:
        for i in range(preruns):
            logger.info(f'Prerun {i} out of {preruns}')
            accept_count = self.sampler.run(preN)
            accept_rate  = accept_count / preN * 100
            logger.info(f'Prerun {i}: acceptance rate is {accept_rate:3.0f}%')
            self.sampler.adapt()
        self.sampler.clear()

    def sample(self, N: int, stride: int) -> np.ndarray:
        sample_total  = N * stride
        sample_chunk  = sample_total // 100
        sample_chunks = [sample_chunk for i in range(0, 99)]
        sample_chunks.append(sample_total - 99 * sample_chunk)
        accept_count  = 0
        for current_chunk in sample_chunks:
            accept_count = accept_count + self.sampler.run(current_chunk)
        accept_rate  = accept_count / (N * stride) * 100
        
        parameter_samples = self.sampler.samples[:][::stride]
        return parameter_samples
    
    def sample_mcmc(self, N: int, stride: int, preruns: int, preN: int, start_point: list[float] = None, cov_scale: float = 0.1) -> np.ndarray:
        self.init_sampler(start_point, cov_scale)
        self.prerun(preruns, preN)
        return self.sample(N, stride)