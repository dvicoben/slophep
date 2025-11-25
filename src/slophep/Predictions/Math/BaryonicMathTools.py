import numpy as np
from flavio.physics.bdecays.common import lambda_K
from flavio.physics import ckm

def LbToLcEllNu_norm_EOS(q2: float, 
                         mB:float, mM: float, ml: float, 
                         par: dict, qi_qj: str) -> float:
    """Calculate decay rate prefactor

    Parameters
    ----------
    q2 : float
        Di-lepton invariant mass
    mB : float
        B-baryon mass
    mM : float
        Final baryon mass
    ml : float
        Charged lepton mass
    par : dict
        flavio parameter dict
    qi_qj : str
        Quark transition

    Returns
    -------
    float
        Rate pre-factor
    """
    lam = lambda_K(mB**2, mM**2, q2)
    GF = par["GF"]
    if qi_qj == 'bu':
        Vij = ckm.get_ckm(par)[0,2] # V_{ub} for b->u transitions
    if qi_qj == 'bc':
        Vij = ckm.get_ckm(par)[1,2] # V_{cb} for b->c transitions
    if q2 <= ml**2 or lam <= 0:
        return 0
    N = GF*Vij*(1 - ml**2/q2)*np.sqrt(q2 / 3.0 / 128 / (np.pi * mB)**3 * np.sqrt(lam))
    return N    


def LbToLcEllNu_amplitudes(q2: float, 
                           mB: float, mM: float, 
                           mb: float, mlight: float, 
                           ml: float, ff: dict, wc: dict, 
                           N: float) -> dict:
    """Calculate decay amplitudes.
    Replicates EOS implementation https://github.com/eos/eos/blob/v1.0.13/eos/b-decays/lambdab-to-lambdac-l-nu.cc.

    Parameters
    ----------
    q2 : float
        Di-lepton invariant mass
    mB : float
        B-baryon mass
    mM : float
        Final baryon mass
    mb : float
        b-quark mass
    mlight : float
        resulting quark mass
    ml : float
        Charged lepton mass
    ff : dict
        Form factors
    wc : dict
        Wilson coefficients
    N : float
        Normalisation prefactor

    Returns
    -------
    dict
        Amplitudes
    """
    beta = 1. - (ml**2/q2)
    m_l_hat = np.sqrt(1 - beta)
    sqrtsminus = np.sqrt((mB - mM)**2 - q2)
    sqrtsplus = np.sqrt((mB + mM)**2 - q2)
    sqrts = np.sqrt(q2)

    # These from using flavio.physics.bdecays.wilsoncoefficents get_wceff_fccc_std
    cvl = wc["VL"]
    cvr = wc["VR"]
    csl = wc["SL"]
    csr = wc["SR"]
    ct = wc["T"]

    # [BKTvD2019] re-implementation of EOS calculation
    A = {
        "perp_1_L" : -2.0 * N * ff["Vp"] * (cvl + cvr) * sqrtsminus,
        "para_1_L" : +2.0 * N * ff["Ap"] * (cvl - cvr) * sqrtsplus,
        "perp_0_L" : +np.sqrt(2.0) * N * ff["V0"] * ((mB + mM) / sqrts) * (cvl + cvr) * sqrtsminus,
        "para_0_L" : -np.sqrt(2.0) * N * ff["A0"] * ((mB - mM) / sqrts) * (cvl - cvr) * sqrtsplus,
        "perp_t_L" : +np.sqrt(2.0) * N * sqrtsplus  * ff["Vt"] * ( 
            ((mB - mM) / sqrts) * (cvl + cvr) 
            + ((mB - mM) / (mb - mlight)) * (csl + csr) / m_l_hat 
        ),
        "para_t_L" : -np.sqrt(2.0) * N * sqrtsminus * ff["At"] * ( 
            ((mB + mM) / sqrts) * (cvl - cvr) 
            - ((mB + mM) / (mb + mlight)) * (csl - csr) / m_l_hat 
        ),
        "para_0_T" : -np.sqrt(8.0) * N * ff["T50"] * sqrtsplus  * ct,
        "perp_0_T" : -np.sqrt(8.0) * N * ff["T0"] * sqrtsminus * ct,
        "para_1_T" : +np.sqrt(4.0) * N * ff["T5p"] * sqrtsplus  * ct * (mB - mM) / sqrts,
        "perp_1_T" : +np.sqrt(4.0) * N * ff["Tp"]  * sqrtsminus * ct * (mB + mM) / sqrts,
    }
    return A


def LbToLcEllNu_observables(q2: float, 
                            mB: float, mM: float, ml: float, 
                            A: dict, alpha: float) -> dict:
    """Calculation of Lb->LcEllNu observables, as per https://arxiv.org/abs/1907.12554.
    Replicates EOS implementation https://github.com/eos/eos/blob/v1.0.13/eos/b-decays/lambdab-to-lambdac-l-nu.cc.

    Parameters
    ----------
    q2 : float
        Di-lepton invariant mass
    mB : float
        B-baryon mass
    mM : float
        Final baryon mass
    ml : float
        Charged lepton mass
    A : dict
        Amplitudes
    alpha : float
        Lc decay parameter

    Returns
    -------
    dict
        Observables
    """
    beta = 1. - (ml**2/q2)
    sqrt1mbeta = np.sqrt(1. - beta)
    K = {}

    K["1ss"] = (
        + 2.0                * (np.abs(A["para_0_L"])**2 + np.abs(A["perp_0_L"])**2)
        + (2.0 - beta)       * (np.abs(A["para_1_L"])**2 + np.abs(A["perp_1_L"])**2 + np.abs(A["para_1_T"])**2 + np.abs(A["perp_1_T"])**2)
        + 2.0 * (1.0 - beta) * (np.abs(A["para_t_L"])**2 + np.abs(A["perp_t_L"])**2 + np.abs(A["para_0_T"])**2 + np.abs(A["perp_0_T"])**2)
        - 4.0 * sqrt1mbeta   * np.real(
            A["para_0_T"] * np.conj(A["para_0_L"]) + A["perp_0_T"] * np.conj(A["perp_0_L"])
            + A["para_1_T"] * np.conj(A["para_1_L"]) + A["perp_1_T"] * np.conj(A["perp_1_L"])
        )
    ) / 4.0
    K["1cc"] = (
        +                    (np.abs(A["para_1_L"])**2 + np.abs(A["perp_1_L"])**2 + np.abs(A["para_0_T"])**2 + np.abs(A["perp_0_T"])**2)
        + (1.0 - beta)     * (np.abs(A["para_0_L"])**2 + np.abs(A["perp_0_L"])**2 + np.abs(A["para_t_L"])**2 + np.abs(A["perp_t_L"])**2 + np.abs(A["para_1_T"])**2 + np.abs(A["perp_1_T"])**2)
        - 2.0 * sqrt1mbeta * np.real(
                A["para_0_T"] * np.conj(A["para_0_L"]) + A["perp_0_T"] * np.conj(A["perp_0_L"])
            + A["para_1_T"] * np.conj(A["para_1_L"]) + A["perp_1_T"] * np.conj(A["perp_1_L"])
        )
    ) / 2.0
    K["1c"] = np.real(
        +                 A["perp_1_L"] * np.conj(A["para_1_L"])
        + (1.0 - beta) * (A["para_0_L"] * np.conj(A["para_t_L"]) + A["perp_0_L"] * np.conj(A["perp_t_L"]) + A["perp_1_T"] * np.conj(A["para_1_T"]))
        - sqrt1mbeta   * (
              A["perp_1_T"] * np.conj(A["para_1_L"]) + A["para_1_T"] * np.conj(A["perp_1_L"])
            + A["para_0_T"] * np.conj(A["para_t_L"]) + A["perp_0_T"] * np.conj(A["perp_t_L"])
        )
    )
    K["2ss"] = alpha * np.real(
        + 2.0                *  A["perp_0_L"] * np.conj(A["para_0_L"])
        + (2.0 - beta)       * (A["perp_1_L"] * np.conj(A["para_1_L"]) - A["perp_1_T"] * np.conj(A["para_1_T"]))
        + 2.0 * (1.0 - beta) * (A["perp_t_L"] * np.conj(A["para_t_L"]) - A["perp_0_T"] * np.conj(A["para_0_T"]))
        - 2.0 * sqrt1mbeta   * (
              A["perp_0_T"] * np.conj(A["para_0_L"]) + A["para_0_T"] * np.conj(A["perp_0_L"])
            + A["perp_1_T"] * np.conj(A["para_1_L"]) + A["para_1_T"] * np.conj(A["perp_1_L"])
        )
    ) / 2.0
    K["2cc"] = alpha * np.real(
        +                (A["perp_1_L"] * np.conj(A["para_1_L"]) + A["perp_0_T"] * np.conj(A["para_0_T"]))
        + (1.0 - beta) * (A["perp_0_L"] * np.conj(A["para_0_L"]) + A["perp_t_L"] * np.conj(A["para_t_L"]) + A["perp_1_T"] * np.conj(A["para_1_T"]))
        - sqrt1mbeta   * (
              A["perp_0_T"] * np.conj(A["para_0_L"]) + A["para_0_T"] * np.conj(A["perp_0_L"])
            + A["perp_1_T"] * np.conj(A["para_1_L"]) + A["para_1_T"] * np.conj(A["perp_1_L"])
        )
    )
    K["2c"] = alpha * (
        +                      (np.abs(A["para_1_L"])**2 + np.abs(A["perp_1_L"])**2)
        + (1.0 - beta)       * (np.abs(A["para_1_T"])**2 + np.abs(A["perp_1_T"])**2)
        + 2.0 * (1.0 - beta) * np.real(
            A["perp_0_L"] * np.conj(A["para_t_L"]) + A["para_0_L"] * np.conj(A["perp_t_L"])
        )
        - 2.0 * sqrt1mbeta * np.real(
              A["para_1_T"] * np.conj(A["para_1_L"]) + A["perp_1_T"] * np.conj(A["perp_1_L"])
            + A["perp_0_T"] * np.conj(A["para_t_L"]) + A["para_0_T"] * np.conj(A["perp_t_L"])
        )
    ) / 2.0
    K["3sc"] = alpha * beta * np.imag(
        + A["perp_1_L"] * np.conj(A["perp_0_L"]) - A["para_1_L"] * np.conj(A["para_0_L"])
        + A["para_1_T"] * np.conj(A["para_0_T"]) - A["perp_1_T"] * np.conj(A["perp_0_T"])
    ) / np.sqrt(2.0)
    K["3s"] = alpha * np.imag(
        + (A["para_1_L"] * np.conj(A["perp_0_L"]) - A["perp_1_L"] * np.conj(A["para_0_L"]))
        + (1.0 - beta) * (
            + A["para_1_L"] * np.conj(A["para_t_L"]) - A["perp_1_L"] * np.conj(A["perp_t_L"])
            + A["para_1_T"] * np.conj(A["perp_0_T"]) - A["perp_1_T"] * np.conj(A["para_0_T"])
        )
        + sqrt1mbeta   * (
            + A["perp_0_T"] * np.conj(A["para_1_L"]) + A["perp_1_T"] * np.conj(A["para_0_L"]) + A["perp_1_T"] * np.conj(A["perp_t_L"])
            - A["para_0_T"] * np.conj(A["perp_1_L"]) - A["para_1_T"] * np.conj(A["perp_0_L"]) - A["para_1_T"] * np.conj(A["para_t_L"])
        )
    ) / np.sqrt(2.0)
    K["4sc"] =  alpha * beta * np.real(
          A["perp_1_L"] * np.conj(A["para_0_L"]) - A["para_1_L"] * np.conj(A["perp_0_L"])
        + A["perp_0_T"] * np.conj(A["para_1_T"]) - A["perp_1_T"] * np.conj(A["para_0_T"])
    ) / np.sqrt(2.0)
    K["4s"] = alpha * np.real(
        + (A["para_1_L"] * np.conj(A["para_0_L"]) - A["perp_1_L"] * np.conj(A["perp_0_L"]))
        + (1.0 - beta) * (
            + A["para_1_L"] * np.conj(A["perp_t_L"]) - A["perp_1_L"] * np.conj(A["para_t_L"])
            + A["para_1_T"] * np.conj(A["para_0_T"]) - A["perp_1_T"] * np.conj(A["perp_0_T"])
        )
        + sqrt1mbeta   * (
            + A["perp_0_T"] * np.conj(A["perp_1_L"]) + A["perp_1_T"] * np.conj(A["perp_0_L"]) + A["perp_1_T"] * np.conj(A["para_t_L"])
            - A["para_0_T"] * np.conj(A["para_1_L"]) - A["para_1_T"] * np.conj(A["para_0_L"]) - A["para_1_T"] * np.conj(A["perp_t_L"])
        )
    ) / np.sqrt(2.0)
    return K