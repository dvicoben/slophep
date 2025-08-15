import numpy as np

def int_one(vmin: float, vmax: float) -> float:
    """Integral of 1"""
    return vmax - vmin

def int_x(vmin: float, vmax: float) -> float:
    "Integral of x"
    return 0.5*(vmax*vmax - vmin*vmin)

def int_x2(vmin: float, vmax: float) -> float:
    "Integral of x^2"
    return (1./3.)*(vmax**3 - vmin**3)

def int_x3(vmin: float, vmax: float) -> float:
    "Integral of x^3"
    return 0.25*(vmax**4 - vmin**4)

def int_x4(vmin: float, vmax: float) -> float:
    "Integral of x^4"
    return 0.2*(vmax**5 - vmin**5)

def int_1minxsq(vmin: float, vmax: float) -> float:
    "Integral of 1-x^2"
    return int_one(vmin, vmax) - int_x2(vmin, vmax)

def int_1plusxsq(vmin: float, vmax: float) -> float:
    "Integral of 1+x^2"
    return int_one(vmin, vmax) + int_x2(vmin, vmax)

def int_xsqmin1(vmin: float, vmax: float) -> float:
    "Integral of x^2-1"
    return int_x2(vmin, vmax) - int_one(vmin, vmax)

def int_2xsqmin1(vmin: float, vmax: float) -> float:
    "Integral of 2x^2-1"
    return 2.*int_x2(vmin, vmax) - int_one(vmin, vmax)

def int_3xsqmin1(vmin: float, vmax: float) -> float:
    "Integral of 3x^2-1"
    return 3.*int_x2(vmin, vmax) - int_one(vmin, vmax)

def int_3xsqmin1sq(vmin: float, vmax: float) -> float:
    "Integral of (3x^2-1)^2"
    return 9.*int_x4(vmin, vmax) - 6.*int_x2(vmin, vmax), int_one(vmin, vmax)

def int_cosx(vmin: float, vmax: float) -> float:
    "Integral of cosx"
    return np.sin(vmax) - np.sin(vmin)

def int_sinx(vmin: float, vmax: float) -> float:
    "Integral of sinx"
    return -np.cos(vmax) + np.cos(vmin)

def int_cos2x(vmin: float, vmax: float) -> float:
    "Integral of cos(2x)"
    return 0.5*(np.sin(2.*vmax) - np.sin(2.*vmin))

def int_sin2x(vmin: float, vmax: float) -> float:
    "Integral of sin(2x)"
    return -0.5*(np.cos(2.*vmax) - np.cos(2.*vmin))

def int_xsqrt1minx(vmin: float, vmax: float) -> float:
    "Integral of x*sqrt(1-x)"
    a=pow(1-vmin,1.5)*(2+3*vmin)
    b=pow(1-vmax,1.5)*(2+3*vmax)
    return (2./15.) * (  a - b  )

def int_2xsqrt1minxsq(vmin: float, vmax: float) -> float:
    "Integral of 2x*sqrt(1-x^2)"
    a = pow(1. - vmin*vmin,1.5)
    b = pow(1. - vmax*vmax,1.5)
    return (2./3.) * (a-b)

def int_sinxwrtcosx(vmin: float, vmax: float) -> float:
    "Integral of sinx w.r.t. cosx"
    return 0.5*(np.sqrt(1-vmax*vmax)*vmax+np.arcsin(vmax)-(np.sqrt(1-vmin*vmin)*vmin+np.arcsin(vmin)))

def int_xn(vmin: float, vmax: float, n: int) -> float:
    "Integral of x^n"
    return pow(vmax,n+1)/(n+1.0)-pow(vmin, n+1)/(n+1.0)
