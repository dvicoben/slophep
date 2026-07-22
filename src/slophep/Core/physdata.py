from __future__ import annotations
import flavio
from typing import Any
import importlib.resources
import json

import logging
logger = logging.getLogger(__file__)


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