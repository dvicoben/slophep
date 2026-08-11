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

from slophep.Observables.ObservableImpl import (
    BToD0stEllNuPrediction, BToD1EllNuPrediction, BToD1stEllNuPrediction, BToD2stEllNuPrediction
)
from slophep.FormFactors.FormFactorsBToDstst import (
    FormFactorBToD0st, FormFactorBToD1, FormFactorBToD1st, FormFactorBToD2st
)


class BuToD1EllNuPrediction(BToD1EllNuPrediction):
    _name = "ObsBuToD1EllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1
                ):
        super().__init__("B+", "D10", lep, nu, FF)


class BuToD0stEllNuPrediction(BToD0stEllNuPrediction):
    _name = "ObsBuToD0stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1
                ):
        super().__init__("B+", "D0*0", lep, nu, FF)


class BuToD1stEllNuPrediction(BToD1stEllNuPrediction):
    _name = "ObsBuToD1stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1st
                ):
        super().__init__("B+", "D1*0", lep, nu, FF)


class BuToD2stEllNuPrediction(BToD2stEllNuPrediction):
    _name = "ObsBuToD2stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD2st
                ):
        super().__init__("B+", "D2*0", lep, nu, FF)



class BdToD1EllNuPrediction(BToD1EllNuPrediction):
    _name = "ObsBdToD1EllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1
                ):
        super().__init__("B0", "D1+", lep, nu, FF)


class BdToD0stEllNuPrediction(BToD0stEllNuPrediction):
    _name = "ObsBdToD0stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD0st,
                ):
        super().__init__("B0", "D0*+", lep, nu, FF)


class BdToD1stEllNuPrediction(BToD1stEllNuPrediction):
    _name = "ObsBdToD1stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1st
                ):
        super().__init__("B0", "D1*+", lep, nu, FF)


class BdToD2stEllNuPrediction(BToD2stEllNuPrediction):
    _name = "ObsBdToD2stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD2st,
                ):
        super().__init__("B0", "D2*+", lep, nu, FF)



class BsToDs1EllNuPrediction(BToD1EllNuPrediction):
    _name = "ObsBsToDs1EllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1,
                ):
        super().__init__("Bs", "Ds1+", lep, nu, FF)


class BsToDs0stEllNuPrediction(BToD0stEllNuPrediction):
    _name = "ObsBsToDs0stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD0st,
                ):
        super().__init__("Bs", "Ds0*+", lep, nu, FF)


class BsToDs1stEllNuPrediction(BToD1stEllNuPrediction):
    _name = "ObsBsToDs1stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD1st,
                ):
        super().__init__("Bs", "Ds1*+", lep, nu, FF)


class BsToDs2stEllNuPrediction(BToD2stEllNuPrediction):
    _name = "ObsBsToDs1stEllNu"
    def __init__(self, 
                 lep: str, 
                 nu: str,
                 FF: FormFactorBToD2st,
                ):
        super().__init__("Bs", "Ds2*+", lep, nu, FF)
