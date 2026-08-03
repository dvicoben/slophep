"""
Module with utilites for ErrorSampler. Contains Fluctuator types and handling of different FluctTypes.
"""
import numpy as np
import scipy.stats as stats
from enum import StrEnum
from typing import Protocol
import scipy.stats as stats

import logging
logger = logging.getLogger(__name__)


class FluctType(StrEnum):
    ARRAYNUMERIC = "arraynumeric",
    DICTNUMERIC  = "dictnumeric",
    NUMERIC      = "numeric"
    NONE         = "None"


def fluctsettings(flucttype: FluctType):
    """Decorator for fluctuatable methods/functions, attaches definition for output-type

    Parameters
    ----------
    flucttype : FluctType
        Function output-type for interaction with ErrorSampler
    """
    def decorate(f):
        setattr(f, "_fluct_type", flucttype)
        return f
    return decorate


def get_sample_err_dictnumeric(samples: list[dict], alpha: float) -> dict[str, tuple[float, float]]:
    merged_res = {k: np.sort([d.get(k) for d in samples]) for k in samples[0]} #set().union(*res)
    # Get upper and lower errors
    errs = {
        k : (ires[int(len(ires)*0.5*alpha)], ires[int(len(ires)*(1-0.5*alpha))]) 
        for k, ires in merged_res.items()
    }
    return errs

def get_sample_err_numeric(samples: list[float], alpha: float) -> tuple[float, float]:
    samples = np.sort(samples)
    return (samples[int(len(samples)*0.5*alpha)], samples[int(len(samples)*(1-0.5*alpha))])

def get_sample_err_arraynumeric(samples: list[list[float]], alpha: float) -> list[tuple[float, float]]:
    samples_tr = np.array(samples).T
    errs = []
    for ientry in samples_tr:
        ierr = get_sample_err_numeric(ientry)
        errs.append(ierr)
    return errs

# Look-up dictionary errorband handling for different output-types as declared in FluctType
_FLUCTTYPE_ERRORCALC = {
    FluctType.NUMERIC      : get_sample_err_numeric,
    FluctType.DICTNUMERIC  : get_sample_err_dictnumeric,
    FluctType.ARRAYNUMERIC : get_sample_err_arraynumeric
}


def check_symmetric(a: np.ndarray, tol: float = 1e-9) -> bool:
    return np.all(np.abs(a-a.T) < tol)

def cov_from_corr(corr: np.ndarray, err_vec: np.ndarray) -> np.ndarray:
    return np.array(corr)*err_vec[:, None]*err_vec[None, :]

def corr_from_cov(cov: np.ndarray) -> np.ndarray:
    cov = np.array(cov)
    err_vec = np.sqrt(np.diag(cov))
    return cov / err_vec[:, None] / err_vec[None, :]


def bifurcated_gaussian_sampler(mean : list[float], 
                                errlo: list[float], 
                                errhi: list[float], 
                                correlation: np.ndarray,
                                rng: np.random.Generator,
                                cutoff: int = -1) -> list[float]:
    """Generate a random sample from bifurcated gaussian using rejection method. 
    For asymmetric errors.

    Parameters
    ----------
    mean : list[float]
        Nominal/mean values for the parameters
    errlo : list[float]
        Lower error
    errhi : list[float]
        Upper error
    correlation : np.ndarray
        Correlation matrix
    rng : np.random.Generator
        Numpy generator object
    cutoff : int, optional
        Number of attempts to cutoff at, by default -1. Negative means no cutoff.

    Returns
    -------
    list[float]
        Fluctuated parameters

    Raises
    ------
    Exception
        If cutoff is reached, raises exception and stops
    """
    attempts = 0
    errlo = np.array(errlo)
    errhi = np.array(errhi)
    while True:
        trial = rng.uniform(-8*errlo, 8*errhi)
        err_vec = np.where(trial < 0, errlo, errhi)
        # trial = np.array([rng.uniform(-8*ilo, 8*ihi) for ilo, ihi in zip(errlo, errhi)])
        # err_vec = np.array([errlo[k] if trial[k] < 0 else errhi[k] for k in range(len(trial))])
        cov = np.array(correlation)*err_vec[:, None]*err_vec[None, :]
        val = stats.multivariate_normal.pdf(trial, mean, cov)
        p = rng.uniform(0,1)
        if p < val:
            return mean+trial
        
        attempts += 1
        if cutoff > 0 and attempts >= cutoff:
            raise Exception("Reached sampling cutoff for bifurcated gaussian")



class Fluctuator(Protocol):
    def generate(self, N: int, rng: np.random.Generator) -> list:
        logger.warning("Blank fluctuator, not generating anything")
        return []


class FlucutatorGaussian:
    def __init__(self, 
                 mean : list[float], 
                 cov  : np.ndarray):
        self._mean: list[float] = mean
        self._cov = np.array(cov)

    @property
    def mean(self) -> list[float]: return self._mean
    @property
    def cov(self) -> np.ndarray: return self._cov

    def generate(self, N: int, rng: np.random.Generator) -> list:
        return rng.multivariate_normal(self.mean, self.cov, int(N))


class FlucutatorBifurGaussian:
    def __init__(self, 
                 mean  : list[float], 
                 errlo : list[float],
                 errhi : list[float],
                 corr  : np.ndarray):
        self._mean  : list[float] = mean
        self._errlo : list[float] = errlo
        self._errhi : list[float] = errhi
        self._corr  : np.ndarray  = np.array(corr)

    @property
    def mean(self) -> list[float]: return self._mean
    @property
    def errlo(self) -> list[float]: return self._errlo
    @property
    def errhi(self) -> list[float]: return self._errhi
    @property
    def corr(self) -> np.ndarray: return self._corr

    def generate(self, N: int, rng: np.random.Generator) -> list:
        fluct = [bifurcated_gaussian_sampler(self.mean, self.errlo, self.errhi, self.corr, rng) for isample in range(N)]
        return fluct

