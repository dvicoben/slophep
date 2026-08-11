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

from typing import Any
from slophep.Core.Parameter import ParameterUser, ParameterManager
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

class FormFactor(ParameterUser):
    _name = "FFBase"

    def set_ff(self, ffpar: dict[str, Any]) -> None:
        for ipar, ival in ffpar.items():
            self.set_userparam(ipar, ival)

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict[str, float]:
        """Calculate form factors at particular q2. To implement in derived class.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            dictionary with FFs 
        """
        raise Exception("get_ff must be implemented in derived class")