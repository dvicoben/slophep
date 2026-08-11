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
from slophep.FormFactors.FormFactorBase import FormFactor

import logging
logger = logging.getLogger(__file__)

class Observable(ParameterUser):
    _name = "ObsBase"
    def __init__(self, 
                 FF: FormFactor):
        """Base class for predictions to handle WC and FF setters"""
        super().__init__()
        self._FF: FormFactor = FF
        pm = self.pm.merge_into(FF.pm)
        # logger.info(f"Manager of {self.name} will be merged into manager of {self.FF.name}.")
        self.set_parammanager(pm)

    @property
    def FF(self) -> FormFactor: 
        """Form factors"""
        return self._FF

    def set_parammanager(self, manager: ParameterManager) -> None:
        super().set_parammanager(manager)
        self.FF.set_parammanager(manager)

    def set_wc(self, wc_dict: dict, eft: str = 'WET', basis: str = 'flavio') -> None:
        """Set the wilson coefficients

        Parameters
        ----------
        wc_dict : dict
            Dictionary of wilson coefficients
        eft : str, optional
            EFT, by default 'WET'
        basis : str, optional
            WC basis (see https://wcxf.github.io/bases.html), by default 'flavio'
        """
        self.pm.set_wc(wc_dict, eft, basis)

    def set_ff(self, ffparams: dict[str, Any]) -> None:
        """Set form factor parameters. Can use None to leave a parameter unchanged.

        Parameters
        ----------
        ffparams : dict
            Dictionary of form factor parameters, names should match FF.ffpars dictionary in the particular 
            scheme being used
        """
        self._FF.set_ff(ffparams)
    
