"""
Constant weight function implementation.

This weight function assigns equal weight to all fragments regardless
of their mass or multiplicity characteristics.
"""

from .weight_function import WeightFunction
from .registry import register_weight_function


@register_weight_function("CONSTANT")
class ConstantWeight(WeightFunction):
    """
    Weight function that assigns equal weight to all fragments.
    
    This is the simplest weight function, giving every fragment the same
    importance in FPIE calculation regardless of its characteristics.
    Typically used when fragment multiplicities are not available or
    when all fragments should be treated equally.
    
    Examples
    --------
    >>> weight_func = ConstantWeight()
    >>> weight_func.compute(100.0, 1)  # Returns 1.0
    1.0
    >>> weight_func.compute(500.0, 5)  # Also returns 1.0
    1.0
    """
    
    def compute(self, mass: float, mult: int) -> float:
        """
        Compute constant weight for any fragment.
        
        Parameters
        ----------
        mass : float
            Fragment mass in Daltons (ignored)
        mult : int
            Fragment multiplicity (ignored)
            
        Returns
        -------
        float
            Always returns 1.0
        """
        return 1.0