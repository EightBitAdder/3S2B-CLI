"""
Exponential decay weight function implementation.

This weight function gives higher weights to fragments with lower
multiplicities, following an exponential decay pattern.
"""

import numpy as np

from .weight_function import WeightFunction
from .registry import register_weight_function


@register_weight_function("EXPONENTIAL-DECAY")
class ExponentialDecayWeight(WeightFunction):
    """
    Weight function using exponential decay based on multiplicity.
    
    Assigns weights according to the formula: exp(-decay_rate * mult)
    This gives higher importance to fragments requiring fewer bond breaks
    (lower multiplicities), as they are generally more likely to form.
    
    Parameters
    ----------
    decay_rate : float, optional
        Rate of exponential decay (default: 0.01)
        Higher values cause faster decay with increasing multiplicity
        
    Examples
    --------
    >>> weight_func = ExponentialDecayWeight()
    >>> weight_func.compute(100.0, 1)  # Returns ~0.99
    0.9900...
    >>> weight_func.compute(100.0, 2)  # Returns ~0.98
    0.9801...
    
    >>> # Custom decay rate
    >>> fast_decay = ExponentialDecayWeight(decay_rate=0.1)
    >>> fast_decay.compute(100.0, 2)  # Returns ~0.82
    0.8187...
    """
    
    def __init__(self, decay_rate: float = 0.01) -> None:
        """
        Initialize exponential decay weight function.
        
        Parameters
        ----------
        decay_rate : float, optional
            Exponential decay rate parameter (default: 0.01)
            
        Raises
        ------
        ValueError
            If decay_rate is negative
        """
        if decay_rate < 0:
            raise ValueError(f"Decay rate must be non-negative, got {decay_rate}")
        self.decay_rate = decay_rate
        
    def compute(self, mass: float, mult: int) -> float:
        """
        Compute exponentially decaying weight based on multiplicity.
        
        Parameters
        ----------
        mass : float
            Fragment mass in Daltons (ignored in this implementation)
        mult : int
            Fragment multiplicity (number of bond breaks)
            
        Returns
        -------
        float
            Weight value following exp(-decay_rate * mult)
        """
        if mult < 0:
            raise ValueError(f"Multiplicity must be non-negative, got {mult}")
        return float(np.exp(-self.decay_rate * mult))
    
    def __repr__(self) -> str:
        """String representation including decay rate."""
        return f"ExponentialDecayWeight(decay_rate={self.decay_rate})"