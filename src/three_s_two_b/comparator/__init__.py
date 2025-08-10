"""
Comparator module for Fragment Pattern Intensity Explained (FPIE) analysis.

This module provides classes and utilities for comparing mass spectral data
against molecular fragment databases using various weight functions.
"""

from .comparator import Comparator, PlotMetaData
from .mass_list import MassList
from .ms_data import MSData

__all__ = [
    "Comparator",
    "PlotMetaData", 
    "MassList",
    "MSData",
]