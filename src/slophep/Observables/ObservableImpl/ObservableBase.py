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
    
