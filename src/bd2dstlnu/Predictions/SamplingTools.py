import numpy as np
import scipy.stats as stats

def bifurcated_gaussian_sampler(mean : list[float], 
                                errlo: list[float], 
                                errhi: list[float], 
                                correlation: np.ndarray,
                                rng: np.random.Generator,
                                cutoff: int = -1):
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
        if attempts >= cutoff:
            raise Exception("Reached sampling cutoff for bifurcated gaussian")