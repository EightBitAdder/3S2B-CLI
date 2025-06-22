from abc import ABC, abstractmethod
import numpy as np


class WeightFunction(ABC):

    @abstractmethod
    def compute(self, mass: float, mult: int) -> float:

        pass