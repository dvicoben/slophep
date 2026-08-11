"""
B0->Pi Form-factors
"""
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
import slophep.FormFactors.FormFactorsBToP as FFBToP
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBdToPi@BSZ"
    def __init__(self):
        super().__init__("B0", "pi+")

    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "f+_0" : 0.1959500850585235,
            "f+_1" : -0.7763565146997043,
            "f+_2" : -0.19660537219243024,
            "f0_1" : -0.01446803146697799,
            "f0_2" : -0.19529667411047308,
            "fT_0" : 0.17311654822482725,
            "fT_1" : -0.40903716524130396,
            "fT_2" : -0.2580355775883392,
            # internalparams
            "m0": 5.540,
            "mp": 5.325
        }
        return ffpar


@FFregistry.register
class BCL(FFBToP.FFBToP_BCL):
    _name = "FFBdToPi@BCL"
    def __init__(self):
        super().__init__("B0", "pi+")