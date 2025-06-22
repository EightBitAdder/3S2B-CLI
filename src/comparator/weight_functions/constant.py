from .weight_function import WeightFunction
from .register import registerWeightFunction


@registerWeightFunction("CONSTANT")
class ConstantWeight(WeightFunction):

    def compute(self, mass: float, mult: int) -> float:

        return 1.0