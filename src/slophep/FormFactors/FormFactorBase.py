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