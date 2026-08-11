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

from __future__ import annotations
import flavio
from typing import Any
import importlib.resources
import json

import logging
logger = logging.getLogger(__name__)


def load_json_data(pkgpath: str, filename: str) -> dict:
    with importlib.resources.open_text(pkgpath, filename) as file:
        data = json.load(file)
        return data


def get_default_params() -> dict:
    par = flavio.default_parameters.get_central_all()
    # Excited D states, from hammer
    extrapar = {
        "m_D10"   : 2.421 ,
        "m_D1+"   : 2.423 ,
        "m_D1*0"  : 2.427 ,
        "m_D1*+"  : 2.427 ,
        "m_D0*0"  : 2.300 ,
        "m_D0*+"  : 2.349 ,
        "m_D2*0"  : 2.461 ,
        "m_D2*+"  : 2.465 ,
        "m_Ds0*+" : 2.3178,
        "m_Ds1*+" : 2.4595,
        "m_Ds1+"  : 2.5351,
        "m_Ds2*+" : 2.5691
    }
    par.update(extrapar)
    return par



DEFAULT_PARAMS = get_default_params()