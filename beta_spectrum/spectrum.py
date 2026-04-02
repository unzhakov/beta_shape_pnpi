# spectrum.py
import numpy as np
from typing import List
from .base import SpectrumComponent


class BetaSpectrum:
    def __init__(self, components: List[SpectrumComponent]):
        self.components = components

    def __call__(self, W: np.ndarray) -> np.ndarray:
        result = np.ones_like(W)

        for comp in self.components:
            result *= comp(W)

        return result
