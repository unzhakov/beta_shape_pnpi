# components/fermi.py
import logging
import numpy as np
from scipy.special import loggamma

from ..base import SpectrumComponent
from ..constants import ALPHA
from ..utils import nuclear_radius


class FermiFunction(SpectrumComponent):
    def __init__(self, Z: int, A: int, logger: logging.Logger | None = None):
        super().__init__(logger=logger)
        self.Z = Z
        self.A = A
        self.alphaZ = ALPHA * Z
        self.gamma = np.sqrt(1.0 - self.alphaZ**2)
        self.R = nuclear_radius(A)
        if self._logger:
            self._logger.debug(
                "FermiFunction: Z=%d, A=%d, alpha*Z=%.4f, gamma=%.6f, R=%.4f",
                Z,
                A,
                self.alphaZ,
                self.gamma,
                self.R,
            )

    def __call__(self, W: np.ndarray) -> np.ndarray:
        if self._logger:
            self._logger.debug(
                "FermiFunction: evaluating at %d points, W range=[%.4f, %.4f]",
                len(W),
                W.min(),
                W.max(),
            )

        p = np.sqrt(W**2 - 1.0)

        # Avoid division by zero
        p = np.clip(p, 1e-30, None)

        eta = self.alphaZ * W / p

        # --- log components ---

        # log prefactor
        log_prefactor = np.log(2 * (1 + self.gamma))

        # log power term
        log_power = -2 * (1 - self.gamma) * np.log(2 * p * self.R)

        # log gamma ratio
        z = self.gamma + 1j * eta
        log_gamma_complex = loggamma(z)
        log_gamma_abs_sq = 2 * np.real(log_gamma_complex)

        log_gamma_denom = 2 * np.real(loggamma(2 * self.gamma + 1))

        # log exponential
        log_exp = np.pi * eta

        logF = log_prefactor + log_power + log_gamma_abs_sq - log_gamma_denom + log_exp

        result = np.asarray(np.exp(logF), dtype=np.float64)

        if self._logger:
            self._logger.debug(
                "FermiFunction: output range=[%.6e, %.6e]",
                result.min(),
                result.max(),
            )

        return result
