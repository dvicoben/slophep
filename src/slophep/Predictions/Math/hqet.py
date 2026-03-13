from math import sqrt, log, pi
from functools import lru_cache
import scipy as sp

def li2(x: float) -> float:
    r"""Complex Dilogarithm"""
    return sp.special.spence(1-x)

def r(w: float) -> float:
    if w == 1:
        return 1
    return log(w + sqrt(-1 + w**2)) / sqrt(-1 + w**2)

def omega_plus(w: float) -> float:
    return w + sqrt(-1 + w**2)


def omega_minus(w: float) -> float:
    return w - sqrt(-1 + w**2)


@lru_cache(maxsize=32)
def omega(w: float, z: float) -> float:
    if w == 1:
        return -1 + (z + 1) / (z - 1) * log(z)
    return (1 + (w * (2 * li2(1 - z * omega_minus(w))
                 - li2(1 - omega_minus(w)**2) -
                 2 * li2(1 - z * omega_plus(w)) + li2(1 - omega_plus(w)**2))
                 ) / (2. * sqrt(-1 + w**2)) - w * log(z) * r(w))


def CS(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (2.*omega(w, z)*(w - wz)*z + (-1. + z**2)*log(z)-(-1. + w)*(z+1.)*(z+1.)*r(w))/(3.*(w - wz)*z)


def CP_Hammer(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (2.*omega(w, z)*(w - wz)*z + (-1. + z**2)*log(z) - (1. + w)*(z+1.)*(z+1.)*r(w))/(3.*(w - wz)*z)


def CP(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return ((-2 * (-w + wz) * (-1 + z) * z *
             (-1 + z + z * (1 + z) * log(z)) *
             ((-1 + z**2) * log(z) +
              (z * (3 + z**2) +
               w *
                 (-1 + z - (3 + 2 * w) * z**2 +
                  z**3)) * r(w)) +
             4 * (w - wz)**2 * z**2 * omega(w, z)) /
            (2. * (w - wz)**2 * z**2))


def CV1(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return ((12 * (-w + wz) * z -
             (-1 + z**2) * log(z) +
             2 * (1 + w) * (-1 + (-1 + 3 * w) * z - z**2) *
             r(w) + 4 * (w - wz) * z * omega(w, z)) /
            (6. * (w - wz) * z))


def CV2(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (-(z *
              (2 * (-w + wz) * (-1 + z) +
               (3 - 2 * w - (-2 + 4 * w) * z + z**2) *
                  log(z)) +
              (2 - (-1 + 5 * w + 2 * w**2) * z +
                  (2 * w + 4 * w**2) * z**2 -
                  (1 + w) * z**3) * r(w)) /
            (6. * (w - wz)**2 * z**2))


def CV3(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return ((2 * (-w + wz) * (-1 + z) * z +
             (1 + (2 - 4 * w) * z + (3 - 2 * w) * z**2) *
             log(z) +
             (1 + w - (2 * w + 4 * w**2) * z +
              (-1 + 5 * w + 2 * w**2) * z**2 - 2 * z**3)
             * r(w)) / (6. * (w - wz)**2 * z))


def CA1(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return ((12 * (-w + wz) * z -
             (-1 + z**2) * log(z) +
             2 * (-1 + w) * (-1 + (1 + 3 * w) * z - z**2) *
             r(w) + 4 * (w - wz) * z * omega(w, z)) /
            (6. * (w - wz) * z))


def CA2(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (-(z *
              (2 * (-w + wz) * (1 + z) +
               (3 + 2 * w - (2 + 4 * w) * z + z**2) *
                  log(z)) +
              (2 + (-1 - 5 * w + 2 * w**2) * z +
                  (-2 * w + 4 * w**2) * z**2 +
                  (1 - w) * z**3) * r(w)) /
            (6. * (w - wz)**2 * z**2))


def CA3(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return ((2 * (-w + wz) * z * (1 + z) -
             (1 - (2 + 4 * w) * z + (3 + 2 * w) * z**2) *
             log(z) +
             (1 - w + (-2 * w + 4 * w**2) * z +
              (-1 - 5 * w + 2 * w**2) * z**2 + 2 * z**3)
             * r(w)) / (6. * (w - wz)**2 * z))


def CT1(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (((-1 + w) *
             (-1 + (2 + 4 * w) * z - z**2) * r(w) +
             (6 * (-w + wz) * z -
              (-1 + z**2) * log(z)) +
             2 * (w - wz) * z * omega(w, z)) / (3. * (w - wz) * z))


def CT2(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (2 * (z * log(z) + (1 - w * z) * r(w))) / (3. * (w - wz) * z)


def CT3(w: float, z: float) -> float:
    wz = 1 / 2 * (z + 1 / z)
    return (2 * (log(z) + (w - z) * r(w))) / (3. * (w - wz))
