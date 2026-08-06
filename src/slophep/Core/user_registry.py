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
