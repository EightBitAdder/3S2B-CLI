"""
Registry system for weight functions.

This module provides registration and retrieval of weight function classes
by string aliases, enabling dynamic selection via CLI and configuration.
"""

from typing import Type, Dict, Optional
import logging

from .weight_function import WeightFunction

logger = logging.getLogger(__name__)

# Global registry mapping aliases to weight function classes
_WEIGHT_FUNCTION_REGISTRY: Dict[str, Type[WeightFunction]] = {}


def register_weight_function(alias: str):
    """
    Decorator to register a weight function class with an alias.
    
    Parameters
    ----------
    alias : str
        String alias for the weight function (case-insensitive)
        
    Returns
    -------
    callable
        Decorator function
        
    Examples
    --------
    >>> @register_weight_function("CUSTOM")
    ... class CustomWeight(WeightFunction):
    ...     def compute(self, mass: float, mult: int) -> float:
    ...         return 1.0
    """
    def decorator(cls: Type[WeightFunction]) -> Type[WeightFunction]:
        if not issubclass(cls, WeightFunction):
            raise TypeError(f"Class {cls.__name__} must inherit from WeightFunction")
            
        alias_upper = alias.upper()
        if alias_upper in _WEIGHT_FUNCTION_REGISTRY:
            logger.warning(f"Weight function alias '{alias}' already registered, overriding")
            
        _WEIGHT_FUNCTION_REGISTRY[alias_upper] = cls
        logger.debug(f"Registered weight function '{alias}' -> {cls.__name__}")
        return cls
    
    return decorator


def fetch_weight_function(alias: Optional[str]) -> WeightFunction:
    """
    Retrieve and instantiate a weight function by alias.
    
    Parameters
    ----------
    alias : str or None
        String alias for the weight function (case-insensitive)
        
    Returns
    -------
    WeightFunction
        Instantiated weight function
        
    Raises
    ------
    ValueError
        If alias is None or not found in registry
        
    Examples
    --------
    >>> weight_func = fetch_weight_function("CONSTANT")
    >>> weight = weight_func.compute(100.0, 2)
    """
    if alias is None:
        raise ValueError("Weight function alias cannot be None")
        
    alias_upper = alias.upper()
    if alias_upper not in _WEIGHT_FUNCTION_REGISTRY:
        available = sorted(_WEIGHT_FUNCTION_REGISTRY.keys())
        raise ValueError(
            f"Unknown weight function '{alias}'. Available options: {available}"
        )
    
    weight_class = _WEIGHT_FUNCTION_REGISTRY[alias_upper]
    return weight_class()


def get_available_weight_functions() -> list[str]:
    """
    Get list of available weight function aliases.
    
    Returns
    -------
    list[str]
        Sorted list of available aliases
    """
    return sorted(_WEIGHT_FUNCTION_REGISTRY.keys())


def clear_registry() -> None:
    """Clear all registered weight functions (mainly for testing)."""
    global _WEIGHT_FUNCTION_REGISTRY
    _WEIGHT_FUNCTION_REGISTRY.clear()
    logger.debug("Cleared weight function registry")