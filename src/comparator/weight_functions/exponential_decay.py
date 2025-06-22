from .weight_function import WeightFunction
from .register import registerWeightFunction
import numpy as np


@registerWeightFunction("EXPONENTIAL-DECAY")
class ExponentialDecayWeight(WeightFunction):

    def compute(self, mass: float, mult: int) -> float:

        return np.exp(-0.01 * mult)