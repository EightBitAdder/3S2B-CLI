"""
Abstract base class for weight functions in FPIE scoring.

Weight functions determine how much each fragment contributes to the
overall FPIE score based on characteristics like mass and multiplicity.
"""

from abc import ABC, abstractmethod


class WeightFunction(ABC):
    """
    Abstract base class for fragment weight functions.
    
    Weight functions take a fragment's mass and multiplicity and return
    a weight value that determines its contribution to the FPIE score.
    Higher weights mean the fragment has more importance in explaining
    the mass spectral data.
    
    Examples
    --------
    >>> class MyWeight(WeightFunction):
    ...     def compute(self, mass: float, mult: int) -> float:
    ...         return 1.0 / mult  # Prefer lower multiplicities
    """
    
    @abstractmethod
    def compute(self, mass: float, mult: int) -> float:
        """
        Compute weight for a fragment.
        
        Parameters
        ----------
        mass : float
            Fragment mass in Daltons
        mult : int
            Fragment multiplicity (number of bond breaks to create)
            
        Returns
        -------
        float
            Weight value (typically between 0 and 1, but not strictly required)
            
        Notes
        -----
        - Weights should be non-negative
        - Higher weights indicate more important fragments
        - Implementation should be deterministic for same inputs
        """
        pass
    
    def __repr__(self) -> str:
        """String representation of weight function."""
        return f"{self.__class__.__name__}()"