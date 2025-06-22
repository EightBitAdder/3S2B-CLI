from typing import Type, Dict
from .weight_function import WeightFunction
import numpy as np


_WEIGHT_FUNCTION_REGISTRY: Dict[str, Type[WeightFunction]] = {}


def registerWeightFunction(alias: str):

    def decorator(cls: Type[WeightFunction]):

        _WEIGHT_FUNCTION_REGISTRY[alias.upper()] = cls

        return cls

    return decorator


def fetchWeightFunction(alias: str) -> WeightFunction:

    try:

        return _WEIGHT_FUNCTION_REGISTRY[alias.upper()]()

    except KeyError:

        raise ValueError()