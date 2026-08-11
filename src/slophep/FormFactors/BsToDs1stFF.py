"""
Bs->Ds1* Form-factors
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

import slophep.FormFactors.FormFactorsBToDstst as FFBToDstst
from slophep.Core.user_registry import FFregistry

@FFregistry.register
class ISGW2(FFBToDstst.FFBToD1st_ISGW2):
    _name = "FFBsToDs1st@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")

    def define_userparams(self):
        ffpar = {
            "msb"     : 5.2                ,
            "msd"     : 0.55               ,
            "bb2"     : 0.54*0.54          ,
            "mbb"     : 5.38               ,
            "msq"     : 1.82               ,
            "bx2"     : 0.41*0.41          ,
            "mbx"     : (3.0*2.54+2.46)/4.0,
            "mqm"     : 0.1                ,
            "nfp"     : 3.0                ,
            "SmearQ2" : True
        }
        return ffpar


@FFregistry.register
class LLSW(FFBToDstst.FFBToD1st_LLSW):
    _name = "FFBsToDs1st@LLSW"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")


@FFregistry.register
class BLR(FFBToDstst.FFBToD1st_BLR):
    _name = "FFBsToDs1st@BLR"
    def __init__(self):
        super().__init__("Bs", "Ds1*+")