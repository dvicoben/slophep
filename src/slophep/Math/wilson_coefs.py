# Copyright (C) 2026  David Vico Benet

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# SLOP or SLOPHEP employs, translates and/or reimplements utilities from:
# - flavio (https://flav-io.github.io/), which is distributed under the MIT License, 
# and without any warranty, see <https://mit-license.org/>
# - Hammer (https://hammer.physics.lbl.gov/), which is distributed under version 3 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>
# - EOS (https://eoshep.org/), which is distributed under version 2 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>

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