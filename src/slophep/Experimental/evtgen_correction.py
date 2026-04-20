import flavio
import numpy as np
from scipy.interpolate import CubicSpline
from typing import Callable

import slophep.Experimental.evtgen_tools as et


def correction_weight_2body(
        mres: float, 
        obs, 
        normscale: float = 1.0
    ) -> float:
    return normscale*mres*et.decrate_int(mres, obs)


def correction_weight_3body(
        mres: float, 
        obs, 
        L: float, 
        width: float, 
        daughters: list[str], 
        normscale: float = 1.0
    ) -> float:
    m_nom = obs.par[f"m_{obs.M}"]
    ls2, ls2_nr = et.threeBody_LS2(mres, m_nom, L, width, daughters, obs.par)
    rate = et.decrate_int(mres, obs)
    return normscale*mres*rate*(ls2/ls2_nr)


def create_interpolated_weightfunc(
        weight_func: Callable,
        obs, 
        mresmin: float, 
        mresmax: float, 
        weight_func_kwargs: dict = {},
        Npoints: int = 1000
    ) -> Callable[[float | np.ndarray], float | np.ndarray]:
    mres   = np.linspace(mresmin, mresmax, int(Npoints))
    # wfunc  = np.vectorize(lambda m: et.decrate_int(m, obs))
    wfunc  = np.vectorize(lambda m: weight_func(m, obs, **weight_func_kwargs))
    gamma  = obs.Gamma()
    spline = CubicSpline(mres, wfunc(mres))
    def correction_func(mass: float | np.ndarray) -> float | np.ndarray:
        return spline(mass)*(1./gamma)
    return correction_func