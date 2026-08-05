from typing import Any

from slophep.Core.Parameter import ParameterUser

class FFregistry:
    _REGISTRY = {}

    @classmethod
    def register(cls, user: type[ParameterUser]):
        if user.name in cls._REGISTRY:
            raise ValueError(f"{user._name} already in registry")
        else:
            cls._REGISTRY[user._name] = user
        return cls

    @classmethod
    def get(cls, name: str, default: Any = None):
        return cls._REGISTRY.get(name, default)