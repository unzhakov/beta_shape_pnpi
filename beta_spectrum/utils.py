# utils.py
import numpy as np

from .constants import HBAR_C_MEV_FM, ME_MEV, R0_FM


def nuclear_radius(A: float, r0: float = R0_FM) -> float:
    R_fm = r0 * A ** (1 / 3)

    return R_fm / (HBAR_C_MEV_FM / ME_MEV)


def momentum(W, rest_mass: float = 1):
    """
    Helper function to get momentum value from energy.

    Args:
        W np.ndarray: Array of full particle energy values in natural units
        rest_mass float: rest mass of the particle in taural units (m_e = 1 by default for electron)

    Returns:
        np.ndarray: array of particle momentum values.
    """
    return np.sqrt(W**2 - rest_mass**2)


def neutrino_energy(W, W0):
    return W0 - W


def W_to_T(W: float | np.ndarray) -> np.ndarray:
    return (W - 1.0) * ME_MEV


def T_to_W(T_MeV: float | np.ndarray) -> np.ndarray:
    return T_MeV / ME_MEV + 1.0
