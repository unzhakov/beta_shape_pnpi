import numpy as np

from beta_spectrum.constants import ALPHA
from beta_spectrum.base import SpectrumComponent
from beta_spectrum.utils import nuclear_radius


class FiniteSizeL0(SpectrumComponent):
    """
    Finite nuclear size correction L0(W)

    Implements low-Z expansion (Hayen et al.)
    Includes (1 + gamma)/2 prefactor
    """

    def __init__(self, Z: int, A: int):
        self.Z = Z
        self.A = A

        self.R = nuclear_radius(A)

        self.alphaZ = ALPHA * Z
        self.gamma = np.sqrt(1 - self.alphaZ**2)

        # beta-minus decay
        self.sign = -1

    def __call__(self, W: np.ndarray) -> np.ndarray:
        # Terms from expansion
        term1 = self.sign * self.alphaZ * W * self.R
        term2 = (13 / 60) * self.alphaZ**2
        term3 = -0.5 * (self.alphaZ * self.R) / W

        return 1 + term1 + term2 + term3


class ChargeDistributionU(SpectrumComponent):
    def __init__(self, Z: int, A: int):
        self.Z = Z
        self.R = nuclear_radius(A)
        self.alphaZ = ALPHA * Z

    def __call__(self, W: np.ndarray) -> np.ndarray:
        return 1.0 + (1 / 5) * (self.alphaZ * W * self.R) ** 2
