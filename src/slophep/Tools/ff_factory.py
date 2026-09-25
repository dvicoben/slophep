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

from typing import Any, Callable
import copy
from slophep.FormFactors.FormFactorBase import FormFactor

class FormFactorFactory:
    @classmethod
    def create(cls, name: str, 
               params: dict[str, Any], 
               ffcalc: Callable[[FormFactor, float], dict[str, float]], 
               base: type[FormFactor] = FormFactor) -> type[FormFactor]:
        """Creates a new FF class, binding params to .define_userparams method, and ffcalc to .calc_ff method.

        Parameters
        ----------
        name : str
            Name of FF scheme
        params : dict[str, Any]
            Dictionary of FF parameters
        ffcalc : Callable[[FormFactor, float], dict[str, float]]
            Callable implementing FF calculation, in manner analogous to calc_ff.
        base : type[FormFactor], optional
            Base class to use for new FF scheme, by default FormFactor

        Returns
        -------
        type[FormFactor]
            The new FF scheme.
        """
        paramd = copy.deepcopy(params) if params is not None else {}
        class newFF(base):
            _name = name
            def define_userparams(self):
                return paramd
            def calc_ff(self, *args, **kwargs):
                return ffcalc(self, *args, **kwargs)

        return newFF