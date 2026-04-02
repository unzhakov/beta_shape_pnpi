# components/fermi.py
import numpy as np
from scipy.special import loggamma

from ..base import SpectrumComponent
from ..constants import ALPHA
from ..utils import nuclear_radius


class FermiFunction(SpectrumComponent):
    def __init__(self, Z: int, A: int):
        self.Z = Z
        self.A = A
        self.alphaZ = ALPHA * Z
        self.gamma = np.sqrt(1.0 - self.alphaZ**2)
        self.R = nuclear_radius(A)

    def __call__(self, W: np.ndarray) -> np.ndarray:
        p = np.sqrt(W**2 - 1.0)

        # Avoid division by zero
        p = np.where(p == 0, 1e-30, p)

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

        log_gamma_denom = 2 * loggamma(2 * self.gamma + 1).real

        # log exponential
        log_exp = np.pi * eta

        logF = (
            log_prefactor
            + log_power
            + log_gamma_abs_sq
            - log_gamma_denom
            + log_exp
        )

        return np.exp(logF)