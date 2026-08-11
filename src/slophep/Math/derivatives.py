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

import numpy as np
import sys

_DER_EPS = 1.e-4

def derivative_left(f, val: float, eps: float = _DER_EPS) -> float:
    return (4.*f(val-eps)-f(val-2.*eps)-3.*f(val))/(-2.*eps)

def derivative_right(f, val: float, eps: float = _DER_EPS) -> float:
    return (4.*f(val+eps)-f(val+2.*eps)-3.*f(val))/(2.*eps)

def derivative_centre(f, val: float, eps: float = _DER_EPS) -> float:
    return (f(val + eps) - f(val - eps))/(2.*eps)

def derivative(f, val: float, eps: float = _DER_EPS, 
               bound_lo: float = -sys.float_info.max, 
               bound_hi: float = sys.float_info.max):
    if val <= bound_lo:
        return derivative_right(f, bound_lo, eps)
    if val >= bound_hi:
        return derivative_left(f, bound_hi, eps)
    
    return derivative_centre(f, val, eps)