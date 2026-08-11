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
from slophep.FormFactors.FormFactorsBToDstst.FFBToDststBase import (
    FormFactorBToD0st, FormFactorBToD1, FormFactorBToD1st, FormFactorBToD2st
)
from slophep.Tools.errfluct_tools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


class FFBToD1_LLSW(FormFactorBToD1):
    _name = "FFBToD1@LLSW"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "t1"   : 0.71,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.0 ,
            "eta2" : 0.0 ,
            "eta3" : 0.0 ,
            #internalparams
            "mb"  : 4.2,
            "mc"  : 1.4,
            "laB" : 0.4,
            "laP" : 0.8 
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)

        eB   = 0.5 / self.get_userparam("mb")
        eC   = 0.5 / self.get_userparam("mc")
        laB  = self.get_userparam("laB")
        laP  = self.get_userparam("laP")
        t1   = self.get_userparam("t1")
        tp   = self.get_userparam("tp")
        tau1 = self.get_userparam("tau1")
        tau2 = self.get_userparam("tau2")
        eta1 = self.get_userparam("eta1")
        eta2 = self.get_userparam("eta2")
        eta3 = self.get_userparam("eta3")

        LambdaD32 = -laB + laP*w
        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        # LOIWtau = t1 + tp*(w-1)  # Hammer
        LOIWtau = t1*(1. + tp*(w-1)) # basf2 EvtGenLLSWFF

        Fv1 = (1 - w*w - eC*(4*(1 + w)*LambdaD32 - (-1 + w*w)*(2*eta1 + 3*eta3 + 3*tau1 - 3*tau2)) - (-1 + w*w)*eB*Fb)/np.sqrt(6.)
        Fv2 = (-3 - eC*(10*eta1 + 4*(w-1)*eta2 - 5*eta3 + (-1 + 4*w)*tau1 + 5*tau2) - 3*eB*Fb)/np.sqrt(6.)
        Fv3 = (-2 + w + eC*(4*LambdaD32 - 2*(6 + w)*eta1 - 4*(w-1)*eta2 - (-2 + 3*w)*eta3 + (2 + w)*tau1 + (2 + 3*w)*tau2) + (2 + w)*eB*Fb)/np.sqrt(6.)
        Fa = (-1 - w - eC*(4*LambdaD32 - (1 + w)*(2*eta1 + 3*eta3) - 3*(w-1)*(tau1 - tau2)) - (w-1)*eB*Fb)/np.sqrt(6.)
        
        ffs = {
            "fS"  : 0.0,
            "fV1" : LOIWtau*Fv1,
            "fV2" : LOIWtau*Fv2,
            "fV3" : LOIWtau*Fv3,
            "fA"  : LOIWtau*Fa ,
            "fT1" : 0.0,
            "fT2" : 0.0,
            "fT3" : 0.0
        }
        return ffs



class FFBToD2st_LLSW(FormFactorBToD2st):
    _name = "FFBToD2st@LLSW"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "t1"   : 0.71,
            "tp"   : -1.6,
            "tau1" : -0.5,
            "tau2" : 2.9 ,
            "eta1" : 0.0 ,
            "eta2" : 0.0 ,
            "eta3" : 0.0 ,
            #internalparams
            "mb"  : 4.2,
            "mc"  : 1.4,
            "laB" : 0.4,
            "laP" : 0.8 
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)

        eB   = 0.5 / self.get_userparam("mb")
        eC   = 0.5 / self.get_userparam("mc")
        laB  = self.get_userparam("laB")
        laP  = self.get_userparam("laP")
        t1   = self.get_userparam("t1")
        tp   = self.get_userparam("tp")
        tau1 = self.get_userparam("tau1")
        tau2 = self.get_userparam("tau2")
        eta1 = self.get_userparam("eta1")
        eta2 = self.get_userparam("eta2")
        eta3 = self.get_userparam("eta3")

        Fb = laB + laP - tau2 - tau1*(1 + 2*w)
        # LOIWtau = t1 + tp*(w-1)  # Hammer
        LOIWtau = t1*(1. + tp*(w-1)) # basf2 EvtGenLLSWFF

        Ka1 = -1. - w - eC*(-((1 + w)*(2*eta1 - eta3)) + (w-1)*(tau1 - tau2)) - (w-1)*eB*Fb
        Ka2 = -2.*eC*(eta2 + tau1)
        Ka3 = 1. - eC*(2*eta1 - 2*eta2 - eta3 + tau1 + tau2) + eB*Fb
        Kv  = -1. - eC*(-2*eta1 + eta3 + tau1 - tau2) - eB*Fb

        ffs = {
            "kP"  : 0.0,
            "kA1" : LOIWtau*Ka1,
            "kA2" : LOIWtau*Ka2,
            "kA3" : LOIWtau*Ka3,
            "kV"  : LOIWtau*Kv ,
            "kT1" : 0.0,
            "kT2" : 0.0,
            "kT3" : 0.0
        }
        return ffs



class FFBToD0st_LLSW(FormFactorBToD0st):
    _name = "FFBToD0st@LLSW"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "zt1"  : 0.68 ,
            "ztp"  : -0.2 ,
            "zeta1": 0.3  ,
            "chi1" : 0.   ,
            "chi2" : 0.   ,
            #internalparams
            "mb"   : 4.2  ,
            "mc"   : 1.4  ,
            "laB"  : 0.4  ,
            "laS"  : 0.76 ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)
        eB = 0.5 / self.get_userparam("mb")
        eC = 0.5 / self.get_userparam("mc")
        zt1 = self.get_userparam("zt1")
        ztp = self.get_userparam("ztp")
        zeta1 = self.get_userparam("zeta1")
        chi1 = self.get_userparam("chi1")
        chi2 = self.get_userparam("chi2")
        laB = self.get_userparam("laB")
        laS = self.get_userparam("laS")

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        # The Belle II EvtGen has callculation of LO IW that differs from the EvtGen code
        # in https://github.com/belle2/basf2/blob/main/generators/evtgen/models/src/EvtLLSWFF.cc
        # LOIWzeta = zt1 + (w-1)*ztp
        LOIWzeta = zt1*(1. + (w-1)*ztp)

        gp = -(eC*((3*LambdaD12)/(1 + w) - 2*(w-1)*zeta1)) - eB*Gb
        gm = 1 + eC*(6*chi1 - 2*(1 + w)*chi2)

        ffs = {
            "gP" : 0.0,
            "g+" : LOIWzeta*gp,
            "g-" : LOIWzeta*gm,
            "gT" : 0.0
        }
        return ffs



class FFBToD1st_LLSW(FormFactorBToD1st):
    _name = "FFBToD1st@LLSW"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "zt1"  : 0.68 ,
            "ztp"  : -0.2 ,
            "zeta1": 0.3  ,
            "chi1" : 0.   ,
            "chi2" : 0.   ,
            #internalparams
            "mb"   : 4.2  ,
            "mc"   : 1.4  ,
            "laB"  : 0.4  ,
            "laS"  : 0.76 ,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        Mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        Mc = mC
        Mb2 = Mb*Mb
        Mc2 = Mc*Mc
        w = (Mb2 + Mc2 - q2)/(2.*Mb*Mc)
        eB = 0.5 / self.get_userparam("mb")
        eC = 0.5 / self.get_userparam("mc")
        zt1 = self.get_userparam("zt1")
        ztp = self.get_userparam("ztp")
        zeta1 = self.get_userparam("zeta1")
        chi1 = self.get_userparam("chi1")
        chi2 = self.get_userparam("chi2")
        laB = self.get_userparam("laB")
        laS = self.get_userparam("laS")

        LambdaD12 = -laB + laS*w
        Gb = (-(laB*(2 + w)) + laS*(1 + 2*w))/(1 + w) - 2*(w-1)*zeta1
        # The Belle II EvtGen has callculation of LO IW that differs from the EvtGen code
        # in https://github.com/belle2/basf2/blob/main/generators/evtgen/models/src/EvtLLSWFF.cc
        # LOIWzeta = zt1 + (w-1)*ztp
        LOIWzeta = zt1*(1. + (w-1)*ztp)

        Gv1 = -1 + w + eC*(LambdaD12 - 2*(w-1)*chi1) - (1 + w)*eB*Gb
        Gv2 = eC*(2*zeta1 - 2*chi2)
        Gv3 = -1 - eC*(LambdaD12/(1 + w) + 2*zeta1 - 2*chi1 + 2*chi2) + eB*Gb
        Ga = 1 + eC*(LambdaD12/(1 + w) - 2*chi1) - eB*Gb

        ffs = {
            "gS"  : 0.0,
            "gV1" : LOIWzeta*Gv1,
            "gV2" : LOIWzeta*Gv2,
            "gV3" : LOIWzeta*Gv3,
            "gA"  : LOIWzeta*Ga,
            "gT1" : 0.0,
            "gT2" : 0.0,
            "gT3" : 0.0
        }

        return ffs