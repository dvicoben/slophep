from math import sqrt
import numpy as np

from slophep.Predictions.Math import integrals as aint

def h_to_f(mB, mP, h, q2):
    """Convert HQET form factors to the standard basis.

    See e.g. arXiv:1309.0301, eq. (31)"""
    ff = {}
    r = mP / mB
    ff['f+'] = ((r + 1) * h['h+'] + (r - 1) * h['h-']) / (2 * sqrt(r))
    fminus = ((r - 1) * h['h+'] + (r + 1) * h['h-']) / (2 * sqrt(r))
    ff['f0'] = ff['f+'] + fminus * q2 / (mB**2 - mP**2)
    ff['fT'] = (r + 1) / (2 * sqrt(r)) * h['hT']
    return ff


def angularPDF(ctl: float, j: dict) -> float:
    return j["a"] + j["b"]*np.cos(ctl) + j["c"]*(np.cos(ctl)**2)


def angular_integrals(ctl_min: float, ctl_max: float) -> dict:
    angint = {
        "a" : aint.int_one(ctl_min, ctl_max),
        "b" : aint.int_x(ctl_min, ctl_max),
        "c" : aint.int_x2(ctl_min, ctl_max)
    }
    return angint

def angularPDF_binned(ctl_min: float, ctl_max: float, j: dict) -> float:
    angint = angular_integrals(ctl_min, ctl_max)
    pdf = j["a"]*angint["a"] + j["b"]*angint["b"] + j["c"]*angint["c"]
    return pdf