import numpy as np

import flavio
from flavio.physics.running import running

def get_CVLSM(par: dict[str, float], 
              scale: float, 
              nf: int):
    r"""Get the Wilson coefficient of the operator $C_V$ in $d_i\to d_j\ell\nu$
    in the SM including EW corrections."""
    if nf >= 4: # for B and D physics
        alpha_e = running.get_alpha(par, scale)['alpha_e']
        return 1 + alpha_e/np.pi * np.log(par['m_Z']/scale)
    else: # for K and pi physics
        # Marciano & Sirlin 1993
        return np.sqrt(1.0232)


def get_wceff_fccc(wc_obj: flavio.WilsonCoefficients, 
                   par: dict[str, float], 
                   qiqj: str, 
                   lep: str, 
                   nu: str, 
                   scale: float, 
                   nf: int = 5,
                   withSM: bool = False):
    r"""Get a dictionary with the $d_i\to d_j$ Wilson coefficients
    in the flavio convention, without or without the SM contribution.
    """
    qqlnu = qiqj + lep + 'nu' + nu
    wc = wc_obj.get_wc(qqlnu, scale, par, nf_out=nf)
    c = {}
    c['VL'] = wc['CVL_'+qqlnu]
    c['VR'] = wc['CVR_'+qqlnu]
    c['SR'] = wc['CSR_'+qqlnu]
    c['SL'] = wc['CSL_'+qqlnu]
    c['T']  = wc['CT_'+qqlnu]
    
    # Whether to add the SM contribution to CVL
    if withSM:
        c_sm = 0.0
        if lep == nu:
            c_sm = get_CVLSM(par, scale, nf)
        c['VL'] += c_sm
    
    return c