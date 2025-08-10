"""
Fragment Pattern Intensity Explained (FPIE) comparator implementation.

This module provides the core comparison logic for calculating FPIE scores
between mass spectral data and molecular fragment databases.
"""

from typing import Union, Optional, Tuple, Dict
from dataclasses import dataclass
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from .weight_functions import WeightFunction, ConstantWeight, ExponentialDecayWeight, fetch_weight_function
from .mass_list import MassList
from .ms_data import MSData

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PlotMetaData:
    """Metadata for FPIE plotting visualization."""
    
    ms_data_masses: np.ndarray
    ms_data_intensities: np.ndarray
    matched_mults: Dict[float, int]
    masses_mults_dict: Dict[float, int]
    matched_weights: Dict[float, float]


class Comparator:
    """
    Compare mass spectral data against molecular fragment lists using FPIE scoring.
    
    The Fragment Pattern Intensity Explained (FPIE) score quantifies how well
    a molecular fragment pattern explains the intensity in mass spectral data.
    
    Parameters
    ----------
    ms_data : MSData
        Mass spectral data containing m/z and intensity values
    mass_list : MassList  
        Molecular fragment mass list with optional multiplicities
    tol : float
        Mass tolerance for matching peaks (in Da)
    weight_function : str, WeightFunction, or None
        Weight function for scoring. If None, defaults based on mass_list.weighted
    
    Examples
    --------
    >>> ms_data = MSData.from_file("spectrum.csv")
    >>> mass_list = MassList.from_file("fragments.csv") 
    >>> comparator = Comparator(ms_data, mass_list, tol=0.1)
    >>> fpie_score, metadata = comparator.calculate_fpie()
    """
    
    def __init__(
        self,
        ms_data: MSData,
        mass_list: MassList,
        *,
        tol: float,
        weight_function: Union[str, WeightFunction, None] = None,
    ) -> None:
        self.ms_data = ms_data
        self.mass_list = mass_list
        self.tol = tol
        self.weight_function = self._resolve_weight_function(weight_function)
        
    @property
    def tol(self) -> float:
        """Mass tolerance in Daltons."""
        return self._tol
    
    @tol.setter
    def tol(self, val: float) -> None:
        if val < 0:
            raise ValueError(f"Tolerance must be non-negative, got {val}")
        self._tol = val
        
    def _resolve_weight_function(self, weight_function: Union[str, WeightFunction, None]) -> WeightFunction:
        """Resolve weight function from various input types."""
        if isinstance(weight_function, WeightFunction):
            return weight_function
        elif isinstance(weight_function, str):
            return fetch_weight_function(weight_function)
        elif weight_function is None:
            # Default based on whether mass list has multiplicities
            return ExponentialDecayWeight() if self.mass_list.weighted else ConstantWeight()
        else:
            raise TypeError(f"Invalid weight function type: {type(weight_function)}")
    
    def _aggregate_ms_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Aggregate MS data by binning masses and summing intensities.
        
        This fixes the critical bug where deduplicated masses were used
        with original intensities array, causing misalignment.
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            Aggregated (masses, intensities) arrays
        """
        # Determine rounding precision based on tolerance
        if self.tol == 0:
            decimals = 0
            rounded_masses = np.round(self.ms_data.masses).astype(int)
        else:
            # Use tolerance to determine appropriate decimal places
            decimals = max(0, int(np.ceil(-np.log10(self.tol))) + 1)
            rounded_masses = np.round(self.ms_data.masses, decimals)
        
        # Create DataFrame for grouping
        df = pd.DataFrame({
            "mz": rounded_masses,
            "intensity": self.ms_data.intensities
        })
        
        # Group by rounded m/z and sum intensities
        aggregated = df.groupby("mz", as_index=False)["intensity"].sum()
        
        return aggregated["mz"].to_numpy(), aggregated["intensity"].to_numpy()