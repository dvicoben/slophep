import flavio
import numpy as np
from scipy.interpolate import CubicSpline
from typing import Callable

def decrate_int(mres: float, obs) -> float:
    mB = obs.par['m_'+obs.B]
    q2min = obs.q2min
    q2max = (mB - mres)**2
    def integrand(qsq):
        rate = obs.dGdq2M(qsq, mres)
        return rate
    return flavio.math.integrate.nintegrate(integrand, q2min, q2max)


def correction_weight(mres: float, obs, normscale: float = 1.0) -> float:
    mB = obs.par['m_'+obs.B]
    q2min = obs.q2min
    q2max = (mB - mres)**2
    def integrand(qsq):
        rate = obs.dGdq2M(qsq, mres)
        return rate
    return normscale*mres*flavio.math.integrate.nintegrate(integrand, q2min, q2max)


def create_interpolated_weightfunc(obs, 
                                   mresmin: float,
                                   mresmax: float,
                                   Npoints: int = 1000
                                   ) -> Callable[[float | np.ndarray], float | np.ndarray]:
    mres   = np.linspace(mresmin, mresmax, int(Npoints))
    wfunc  = np.vectorize(lambda m: decrate_int(m, obs))
    gamma  = obs.Gamma()
    spline = CubicSpline(mres, wfunc(mres))
    def correction_func(mass: float | np.ndarray) -> float | np.ndarray:
        return mass*spline(mass)*(1./gamma)
    return correction_func