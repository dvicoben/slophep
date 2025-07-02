import numpy as np

def calc_unaing_obs(j: dict[float]) -> float:
    norm = 3/4. * (2 * j['1s'] + j['1c']) - 1/4. * (2 * j['2s'] + j['2c'])
    obs = {
        "FL" : (0.25*(3*j["1c"] - j["2c"]))/norm,
        "AFB" : ((3.0/8.0)*(j["6c"] + 2*j["6s"]))/norm,
        "FLt" : 0.25*(j["1c"] - 3*j["2c"] + 2*j["1s"] - 6*j["2s"])/norm,
        3 : j[3],
        9 : j[9]
    }
    return obs

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



def angularPDF(ctx: float, ctl: float, chi: float, j: dict[float]) -> float:
    ctx2 = ctx*ctx
    stx2 = 1 - ctx2
    # tx = np.arccos(ctx)
    # s2tx = np.sin(2*tx)
    stx = np.sqrt(stx2)
    s2tx = 2*stx*ctx

    ctl2 = ctl*ctl
    stl2 = 1 - ctl2
    c2tl = 2*ctl2 - 1
    # tl = np.arccos(ctl)
    # stl = np.sin(tl)
    # s2tl = np.sin(2*tl)
    stl = np.sqrt(stl2)
    s2tl = 2*stl*ctl

    cph = np.cos(chi)
    sph = np.sin(chi)
    c2ph = 2*cph*cph - 1
    s2ph = 2*sph*cph
    
    pdf = 9/(32*np.pi) * (
          j['1c']*ctx2         # I1c * cos2(thetaX)
        + j['1s']*stx2         # I1s * sin2(thetaX)
        + j['2c']*ctx2*c2tl    # I2c * cos2(thetaX) * cos(2thetaL)
        + j['2s']*stx2*c2tl    # I2s * sin2(thetaX) * cos(2thetaL)
        + j[3]*stx2*stl2*c2ph  # I3  * sin2(thetaX) * sin2(thetaL) * cos(2phi)
        + j[4]*s2tx*s2tl*cph   # I4  * sin(2thetaX) * sin(2thetaL) * cos(phi)
        + j[5]*s2tx*stl*cph    # I5  * sin(2thetaX) * sin(thetaL)  * cos(phi)
        + j['6c']*ctx2*ctl     # I6c * cos2(thetaX) * cos(thetaL)
        + j['6s']*stx2*ctl     # I6s * sin2(thetaX) * cos(thetaL)
        + j[7]*s2tx*stl*sph    # I7  * sin(2thetaX) * sin(thetaL)  * sin(phi)
        + j[8]*s2tx*s2tl*sph   # I8  * sin(2thetaX) * sin(2thetaL) * sin(phi)
        + j[9]*stx2*stl2*s2ph  # I9  * sin2(thetaX) * sin2(thetaL) * sin2(phi)
    )
    return pdf


def angular_integrals(
        ctx_min: float, ctx_max: float,
        ctl_min: float, ctl_max: float,
        chi_min: float, chi_max: float) -> dict[float]:
    angint = {
        '1c' : int_x2(ctx_min, ctx_max)*int_one(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        '1s' : int_1minxsq(ctx_min, ctx_max)*int_one(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        '2c' : int_x2(ctx_min, ctx_max)*int_2xsqmin1(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        '2s' : int_1minxsq(ctx_min, ctx_max)*int_2xsqmin1(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        3    : int_1minxsq(ctx_min, ctx_max)*int_1minxsq(ctl_min, ctl_max)*int_cos2x(chi_min, chi_max),
        4    : int_2xsqrt1minxsq(ctx_min, ctx_max)*int_2xsqrt1minxsq(ctl_min, ctl_max)*int_cosx(chi_min, chi_max),
        5    : int_2xsqrt1minxsq(ctx_min, ctx_max)*int_sinxwrtcosx(ctl_min, ctl_max)*int_cosx(chi_min, chi_max),
        '6c' : int_x2(ctx_min, ctx_max)*int_x(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        '6s' : int_1minxsq(ctx_min, ctx_max)*int_x(ctl_min, ctl_max)*int_one(chi_min, chi_max),
        7    : int_2xsqrt1minxsq(ctx_min, ctx_max)*int_sinxwrtcosx(ctl_min, ctl_max)*int_sinx(chi_min, chi_max),
        8    : int_2xsqrt1minxsq(ctx_min, ctx_max)*int_2xsqrt1minxsq(ctl_min, ctl_max)*int_sinx(chi_min, chi_max),
        9    : int_1minxsq(ctx_min, ctx_max)*int_1minxsq(ctl_min, ctl_max)*int_sin2x(chi_min, chi_max)
    }
    return angint


def angularPDF_binned(
        ctx_min: float, ctx_max: float,
        ctl_min: float, ctl_max: float,
        chi_min: float, chi_max: float,
        j: dict[float]) -> float:
    angint = angular_integrals(ctx_min, ctx_max, ctl_min, ctl_max, chi_min, chi_max)
    pdf = 9/(32*np.pi) * (
              j['1c'] * angint['1c']
            + j['1s'] * angint['1s']
            + j['2c'] * angint['2c']
            + j['2s'] * angint['2s']
            + j[3]    * angint[3]
            + j[4]    * angint[4]
            + j[5]    * angint[5]
            + j['6c'] * angint['6c']
            + j['6s'] * angint['6s']
            + j[7]    * angint[7]
            + j[8]    * angint[8]
            + j[9]    * angint[9]
        )
    return pdf

