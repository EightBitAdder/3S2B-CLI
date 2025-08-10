"""
Weight functions for Fragment Pattern Intensity Explained (FPIE) scoring.

This module provides various weight functions that can be applied during
FPIE calculation to adjust the importance of fragments based on their
mass and multiplicity characteristics.
"""

from .weight_function import WeightFunction
from .constant import ConstantWeight
from .exponential_decay import ExponentialDecayWeight
from .registry import register_weight_function, fetch_weight_function

__all__ = [
    "WeightFunction",
    "ConstantWeight", 
    "ExponentialDecayWeight",
    "register_weight_function",
    "fetch_weight_function",
]