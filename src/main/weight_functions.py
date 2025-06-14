from typing import Dict
from abc import ABC, abstractmethod
import numpy as np


_WEIGHT_FUNCTION_REGISTRY: dict[str, type["WeightFunction"]] = {}


def register_weight_function(alias: str):

    def decorator(cls: type["WeightFunction"]):

        _WEIGHT_FUNCTION_REGISTRY[alias.upper()] = cls

        return cls

    return decorator


class WeightFunction(ABC):

    @abstractmethod
    def compute(self, mass: float, mult: int) -> float:

        pass


@register_weight_function("CONSTANT")
class ConstantWeight(WeightFunction):

    def compute(self, mass: float, mult: int) -> float:

        return 1.0


@register_weight_function("EXPONENTIAL-DECAY")
class ExponentialDecayWeight(WeightFunction):

    def compute(self, mass: float, mult: int) -> float:

        return np.exp(-0.01 * mult)


def lsWeightFunctions() -> list[str]:

    return list(_WEIGHT_FUNCTION_REGISTRY.keys())


def fetchWeightFunction(alias: str) -> WeightFunction:

    try:

        return _WEIGHT_FUNCTION_REGISTRY[alias.upper()]()

    except KeyError:

        raise ValueError()