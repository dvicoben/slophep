"""
B+->D0 Form-factors
"""
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
import slophep.FormFactors.FormFactorsBToP as FFBToP
from slophep.Tools.errfluct_tools import FluctType, fluctsettings
from slophep.Core.user_registry import FFregistry

import logging
logger = logging.getLogger(__name__)

@FFregistry.register
class BSZ(FFBToP.FFBToP_BSZ):
    _name = "FFBuToD@BSZ"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BLPR(FFBToP.FFBToP_BLPR):
    _name = "FFBuToD@BLPR"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BLPRXP(FFBToP.FFBToP_BLPRXP):
    _name = "FFBuToD@BLPRXP"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class BGL(FFBToP.FFBToP_BGL):
    _name = "FFBuToD@BGL"
    def __init__(self):
        super().__init__("B+", "D0")


@FFregistry.register
class CLN(FFBToP.FFBToP_CLN):
    _name = "FFBuToD@CLN"
    def __init__(self):
        super().__init__("B+", "D0")


# Hammer equivalent classes
@FFregistry.register
class CLN_Hammer(FFBToP.FFBToP_CLN):
    _name = "FFBuToD@CLN_Hammer"
    def __init__(self):
        super().__init__("B+", "D0")

    def define_userparams(self):
        ffpar = {
            "RhoSq"  : 1.186,
            "G1"     : 1.082,
            "Delta"  : 1.0  ,
            "a"      : 1.0  ,  # zero recoil expansion
            "mcOnmb" : 0.29 ,
            "lamBar" : 0.48 ,
            "ash"    : 0.26/np.pi
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float):
        """B->P Hammer basis differs in fT from standard SLOP basis, use with care for predictions."""
        scale = self.get_param(f"m_{self.B}") + self.get_param(f"m_{self.P}")
        ffs = super().get_ff(q2)
        ffs["fT"] *= 1./scale
        return ffs


@FFregistry.register
class BGL_Hammer(FFBToP.FFBToP_BGLGeneric):
    _name = "FFBuToD@BGL_Hammer"
    def __init__(self):
        super().__init__("B+", "D0", 4, 4)

    def define_userparams(self):
        ffpar = super().define_userparams()
        ffpar.update({
            "a_f+_0" : 0.01565,
            "a_f+_1" : -0.0353,
            "a_f+_2" : -0.043,
            "a_f+_3" : 0.194,
            "a_f0_0" : 0.07932,
            "a_f0_1" : -0.214,
            "a_f0_2" : 0.17,
            "a_f0_3" : -0.958,
            #internalparams
            "BcStatesp" : np.array([6.329, 6.920, 7.020]),
            "BcStates0" : np.array([6.716, 7.121]),
            "ChiT"      : 5.131e-4, # 1606.08030
            "ChiL"      : 6.332e-3, # 1606.08030
            "nmax"      : 4,
            "nc"        : 2.6
        })
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float):
        """B->P Hammer basis differs in fT from standard SLOP basis, use with care for predictions - SM only so shouldn't matter"""
        scale = self.get_param(f"m_{self.B}") + self.get_param(f"m_{self.P}")
        ffs = super().get_ff(q2)
        ffs["fT"] *= 1./scale
        return ffs
    

@FFregistry.register
class BLPR_Hammer(FFBToP.FFBToP_BLPR):
    _name = "FFBuToD@BLPR_Hammer"
    def __init__(self):
        logger.info("B->P BLPR_Hammer basis differs in fT from standard SLOP basis, use with care for predictions.")
        super().__init__("B+", "D0")

    def define_userparams(self):
        ffpar = {
            "a"         : 1.509/np.sqrt(2),
            "rD"        : 1867./5280.,
            "RhoSq"     : 1.24,
            "Chi21"     : -0.06,
            "Chi2p"     : 0.0,
            "Chi3p"     : 0.05,
            "Eta1"      : 0.30,
            "Etap"      : -0.05,
            "dV20"      : 0.0,
            "mb"        : 4.710,
            "delta_mbc" : 3.4,
            "normscale" : 1.0,
            "ash"       : 0.26/np.pi,
            "mbarB"     : 5.313,
            "lam1"      : -0.3 , # GeV^2
            "ebReb"     : 0.861,
            "ecRec"     : 0.822,
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float):
        """B->P Hammer basis differs in fT from standard SLOP basis, use with care for predictions."""
        scale = self.get_param(f"m_{self.B}") + self.get_param(f"m_{self.P}")
        ffs = super().get_ff(q2)
        ffs["fT"] *= 1./scale
        return ffs