# Copyright (C) 2026  David Vico Benet

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# SLOP or SLOPHEP employs, translates and/or reimplements utilities from:
# - flavio (https://flav-io.github.io/), which is distributed under the MIT License, 
# and without any warranty, see <https://mit-license.org/>
# - Hammer (https://hammer.physics.lbl.gov/), which is distributed under version 3 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>
# - EOS (https://eoshep.org/), which is distributed under version 2 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>

from typing import Callable
import numpy as np
import pypmc

import logging
logger = logging.getLogger(__name__)

SamplingFunction = Callable[..., float]

class MCGenerator:
    def __init__(self, seed: int = None):
        """MC Generator

        Parameters
        ----------
        seed : int, optional
            _description_, by default None
        """
        self._seed    : int                          = seed
        self._pdf     : SamplingFunction             = None
        self._logpdf  : SamplingFunction             = None
        self._bounds  : list[tuple[float, float]]    = []
        self._rng     : np.random.RandomState        = np.random.mtrand if seed is None else np.random.mtrand.RandomState(seed)
        self._sampler : pypmc.sampler.markov_chain.AdaptiveMarkovChain = None

    @property
    def pdf(self) -> SamplingFunction:
        """PDF to sample"""
        return self._pdf
    @property
    def bounds(self) -> list[tuple[float, float]]:
        """Bounds of pdf arguments"""
        return self._bounds
    @property
    def seed(self) -> int:
        """rng seed"""
        return self._seed
    @property
    def rng(self) -> np.random.RandomState:
        """Random state object"""
        return self._rng
    @property
    def sampler(self) -> pypmc.sampler.markov_chain.AdaptiveMarkovChain:
        """MCMC sampler"""
        return self._sampler

    def logpdf(self, x: list[float]) -> float:
        if self._logpdf is not None:
            return self._logpdf(*x)
        return np.log(self.pdf(*x))

    def set_logPDF(self, logpdf: SamplingFunction, bounds: list[tuple[float, float]]) -> None:
        """Set the pdf to sample

        Parameters
        ----------
        logpdf : SamplingFunction
            logPDF callable to sample, should have float arguments logpdf(x: float, y: float, z: float)
        bounds : list[tuple[float, float]]
            Bounds for arguments of pdf, in positional order
        """
        self._logpdf = logpdf
        self._bounds = bounds
    
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
        """Prepare the sampler

        Parameters
        ----------
        start_point : list[float], optional
            Sampler starting point, by default None. If None uses random point uniformly sampled from bounds.
        cov_scale : float, optional
            Scaling of proposal cov, by default 0.1
        """
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
        """Pre-run/burnin of sampler

        Parameters
        ----------
        preruns : int
            Number of preruns
        preN : int
            Number of events per prerun
        """
        for i in range(preruns):
            logger.info(f'Prerun {i} out of {preruns}')
            accept_count = self.sampler.run(preN)
            accept_rate  = accept_count / preN * 100
            logger.info(f'Prerun {i}: acceptance rate is {accept_rate:3.0f}%')
            self.sampler.adapt()
        self.sampler.clear()

    def sample(self, N: int, stride: int) -> np.ndarray:
        """Generate samples

        Parameters
        ----------
        N : int
            Number of samples to return
        stride : int
            Number by which the actual amount of samples is thinned to return N samples

        Returns
        -------
        np.ndarray
            The pseudo-data
        """
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
        """Executes init_sampler, prerun, and sample.

        Parameters
        ----------
        N : int
            Number of samples to return
        stride : int
            Number by which the actual amount of samples is thinned to return N samples
        preruns : int
            Number of preruns
        preN : int
            Number of events per prerun
        start_point : list[float], optional
            Sampler starting point, by default None. If None uses random point uniformly sampled from bounds.
        cov_scale : float, optional
            Scaling of proposal cov, by default 0.1

        Returns
        -------
        np.ndarray
            The pseudo-data
        """
        self.init_sampler(start_point, cov_scale)
        self.prerun(preruns, preN)
        return self.sample(N, stride)