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

from math import sqrt
import numpy as np

from slophep.Math import integrals as aint

def h_to_f(mB: float, mP: float, h: dict[str, float], q2: float):
    """Convert HQET form factors to the standard basis.

    See e.g. arXiv:1309.0301, eq. (31)"""
    ff = {}
    r = mP / mB
    ff['f+'] = ((r + 1) * h['h+'] + (r - 1) * h['h-']) / (2 * sqrt(r))
    fminus = ((r - 1) * h['h+'] + (r + 1) * h['h-']) / (2 * sqrt(r))
    ff['f0'] = ff['f+'] + fminus * q2 / (mB**2 - mP**2)
    ff['fT'] = (r + 1) / (2 * sqrt(r)) * h['hT']
    return ff


def angularPDF(ctl: float, j: dict[str, float]) -> float:
    return j["a"] + j["b"]*np.cos(ctl) + j["c"]*(np.cos(ctl)**2)


def angular_integrals(ctl_min: float, ctl_max: float) -> dict[str, float]:
    angint = {
        "a" : aint.int_one(ctl_min, ctl_max),
        "b" : aint.int_x(ctl_min, ctl_max),
        "c" : aint.int_x2(ctl_min, ctl_max)
    }
    return angint

def angularPDF_binned(ctl_min: float, ctl_max: float, j: dict[str, float]) -> float:
    angint = angular_integrals(ctl_min, ctl_max)
    pdf = j["a"]*angint["a"] + j["b"]*angint["b"] + j["c"]*angint["c"]
    return pdf