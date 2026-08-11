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

from typing import Any

from slophep.Core.Parameter import ParameterUser

class UserRegistry:
    def __init__(self):
        self._REGISTRY: dict[str, type[ParameterUser]] = {}

    @property
    def REGISTRY(self) -> dict[str, type[ParameterUser]]: return self._REGISTRY

    def register(self, user: type[ParameterUser]):
        if user.name in self.REGISTRY:
            raise ValueError(f"{user._name} already in registry")
        else:
            self.REGISTRY[user._name] = user
        return user

    def get(self, name: str, default: Any = None):
        return self.REGISTRY.get(name, default)

    def get_instance_registry(self, *args, **kwargs) -> dict[str, ParameterUser]:
        return {iname : iuser(*args, **kwargs) for iname, iuser in self.REGISTRY.items()}

    def find_first_containing(self, snippet: str) -> type[ParameterUser] | None:
        for iname, icls in self.REGISTRY.items():
            if snippet in iname:
                return icls
        return None

    def find_all_containing(self, snippet: str) -> dict[str, type[ParameterUser]]:
        users = {}
        for iname, icls in self.REGISTRY.items():
            if snippet in iname:
                users[iname] = icls
        return users


FFregistry = UserRegistry()
