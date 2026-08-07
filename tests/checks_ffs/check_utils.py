import numpy as np
from os.path import abspath, dirname

from slophep.FormFactors.FormFactorBase import FormFactor

_RELATIVE_TOLERANCE = 1e-4

def within_tolerance(val: np.ndarray, ref: np.ndarray, tol: float = _RELATIVE_TOLERANCE) -> bool:
    return np.allclose(val, ref, rtol=tol)

# Getting the SLOP predictions:
def get_spectrum_slop(ff: FormFactor, 
                      qsq: list[float], 
                      method: str):
    res = {}
    for iq2 in qsq:
        ffmethod = getattr(ff, method)
        iff = ffmethod(iq2)
        for ielem in iff:
            if ielem not in res:
                res[ielem] = []
            res[ielem].append(iff[ielem])
    return {k : np.array(res[k]) for k in res}


def read_txt_data(inpath: str,
                  keyorder: list[str],
                  inpath_rel: bool = True) -> tuple[list[float], dict[str, list[float]]]:
    if inpath_rel:
        script_dir = dirname(abspath(__file__))
        inpath = f"{script_dir}/"+inpath
    data = np.loadtxt(inpath, float, skiprows=1).T
    qsq = data[0]
    spectrum = {keyorder[i] : data[i+1] for i in range(len(keyorder))}
    return qsq, spectrum


def convert_btov_hammer_bglbasis(FFs: dict[str, float], Mb: float, Mc: float, qsq: np.ndarray) -> dict[str, float]:
    Mb2 = Mb**2
    Mb3 = Mb**3
    rC = Mc/Mb
    rC2 = rC**2
    w = (Mb**2 + Mc**2 - qsq) / (2 * Mb * Mc)
    w2 = w**2
    Fpf = (rC-w)/(2.*rC*Mb2*(w2 - 1))
    FpF1 = 1./(2.*rC*Mb3*(w2 - 1))
    Fmf = (rC+w)/(2.*rC*Mb2*(w2 - 1))
    FmF1 = 1./(2.*rC*Mb3*(w2 - 1))*(rC2-1)/(1 + rC2 - 2*rC*w)
    FmP1 = np.sqrt(rC)*(rC+1)/(Mb*(1 + rC2 - 2*rC*w))
    F1 = (1./FpF1)*(FFs["Fp"] - Fpf*FFs["Ff"])
    F2 = (1./FmP1)*((1+rC)/np.sqrt(rC))*(FFs["Fm"] - Fmf*FFs["Ff"] - FmF1*F1)
    ffconv = {
        "F1" : F1,
        "F2" : F2,
        "g"  : FFs["Fg"]*2.,
        "f"  : FFs["Ff"]
    }
    return ffconv


def convert_btov_hammer_hbasis(FFs: dict[str, float], Mb: float, Mc: float, qsq: np.ndarray) -> dict[str, float]:
    ff = {}
    sqMbMc = np.sqrt(Mb*Mc)
    Mb3 = Mb**3
    
    ff["hA1"] = (2. * sqMbMc)*FFs["Ff"]/((Mb + Mc)*(Mb + Mc) - qsq)
    ff["hV"] = (2. * sqMbMc)*FFs["Fg"]
    ff["hA2"] = -(FFs["Fm"] + FFs["Fp"])/np.sqrt(Mc / Mb3)
    ff["hA3"] = 0.5*(2. * sqMbMc)*(FFs["Fm"] - FFs["Fp"])
    ff["hT3"] = FFs["Fzt"]*(2. * pow(Mb * Mc, 1.5))/Mc
    ff["hT2"] = (2. * sqMbMc)*((Mb+Mc)*FFs["Fmt"] + (Mb-Mc)*FFs["Fpt"])/((Mb-Mc)**2 - (Mb+Mc)**2)
    ff["hT1"] = (FFs["Fmt"] + (ff["hT2"] * (Mb + Mc)) / (2. * sqMbMc))*((2. * sqMbMc)/(Mb - Mc))
    return ff

