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
# - MCAmbulance (https://github.com/Herren/MCAmbulance), which is distributed under version 3 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>

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
    m_nom = obs.get_param(f"m_{obs.M}")
    ls2, ls2_nr = et.threeBody_LS2(mres, m_nom, L, width, daughters, obs.pm)
    rate = et.decrate_int(mres, obs)
    return normscale*mres*rate*(ls2/ls2_nr)


def create_interpolated_weightfunc(
        weight_func: Callable,
        obs, 
        mresmin: float, 
        mresmax: float, 
        weight_func_kwargs: dict = None,
        Npoints: int = 1000
    ) -> Callable[[float | np.ndarray], float | np.ndarray]:
    mres   = np.linspace(mresmin, mresmax, int(Npoints))
    # wfunc  = np.vectorize(lambda m: et.decrate_int(m, obs))
    weight_func_kwargs = weight_func_kwargs if weight_func_kwargs is not None else {}
    wfunc  = np.vectorize(lambda m: weight_func(m, obs, **weight_func_kwargs))
    spline = CubicSpline(mres, wfunc(mres))
    def correction_func(mass: float | np.ndarray) -> float | np.ndarray:
        return spline(mass)
    return correction_func