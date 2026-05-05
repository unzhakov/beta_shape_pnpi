import logging

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

    def __init__(self, Z: int, A: int, logger: logging.Logger | None = None):
        super().__init__(logger=logger)
        self.Z = Z
        self.A = A

        self.R = nuclear_radius(A)

        self.alphaZ = ALPHA * Z
        self.gamma = np.sqrt(1 - self.alphaZ**2)

        # beta-minus decay
        self.sign = -1

        if self._logger:
            self._logger.debug(
                "FiniteSizeL0: Z=%d, A=%d, R=%.4f, alphaZ=%.4f, gamma=%.6f",
                Z,
                A,
                self.R,
                self.alphaZ,
                self.gamma,
            )

    def __call__(self, W: np.ndarray) -> np.ndarray:
        if self._logger:
            self._logger.debug(
                "FiniteSizeL0: evaluating at %d points, W range=[%.4f, %.4f]",
                len(W),
                W.min(),
                W.max(),
            )

        # Terms from expansion
        term1 = self.sign * self.alphaZ * W * self.R
        term2 = (13 / 60) * self.alphaZ**2
        term3 = -0.5 * (self.alphaZ * self.R) / W

        result = 1 + term1 + term2 + term3

        if self._logger:
            self._logger.debug(
                "FiniteSizeL0: output range=[%.6e, %.6e]",
                result.min(),
                result.max(),
            )

        return result


class ChargeDistributionU(SpectrumComponent):
    def __init__(self, Z: int, A: int, logger: logging.Logger | None = None):
        super().__init__(logger=logger)
        self.Z = Z
        self.R = nuclear_radius(A)
        self.alphaZ = ALPHA * Z

        if self._logger:
            self._logger.debug(
                "ChargeDistributionU: Z=%d, A=%d, R=%.4f, alphaZ=%.4f",
                Z,
                A,
                self.R,
                self.alphaZ,
            )

    def __call__(self, W: np.ndarray) -> np.ndarray:
        if self._logger:
            self._logger.debug(
                "ChargeDistributionU: evaluating at %d points, W range=[%.4f, %.4f]",
                len(W),
                W.min(),
                W.max(),
            )

        result = 1.0 + (1 / 5) * (self.alphaZ * W * self.R) ** 2

        if self._logger:
            self._logger.debug(
                "ChargeDistributionU: output range=[%.6e, %.6e]",
                result.min(),
                result.max(),
            )

        return result
