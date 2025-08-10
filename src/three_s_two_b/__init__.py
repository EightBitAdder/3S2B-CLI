"""
3S2B CLI and Library for Molecular Fragment Database Analysis.

A command-line interface and Python library implementing the 3S2B algorithm
for Fragment Pattern Intensity Explained (FPIE) scoring and molecular fragment
database analysis.
"""

from typing import TYPE_CHECKING

__version__ = "1.0.0"
__author__ = "Jesse Fraser, Arun Moorthy"
__email__ = "crafts.lab@trentu.ca"

# Public API exports
if TYPE_CHECKING:
    from .comparator import Comparator, MassList, MSData
    from .fragmentor import Fragmentor, fetch_mass_list

__all__ = [
    "__version__",
    "Comparator",
    "MassList", 
    "MSData",
    "Fragmentor",
    "fetch_mass_list",
]


def __getattr__(name: str):
    """Lazy import for public API."""
    if name == "Comparator":
        from .comparator import Comparator
        return Comparator
    elif name == "MassList":
        from .comparator import MassList
        return MassList
    elif name == "MSData":
        from .comparator import MSData
        return MSData
    elif name == "Fragmentor":
        from .fragmentor import Fragmentor
        return Fragmentor
    elif name == "fetch_mass_list":
        from .fragmentor import fetch_mass_list
        return fetch_mass_list
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name__}'")