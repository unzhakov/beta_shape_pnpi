# utils.py
from typing import Union

import numpy as np

from .constants import HBAR_C_MEV_FM, ME_MEV, R0_FM

Real = Union[float, np.ndarray]


def nuclear_radius(A: float, r0: float = R0_FM) -> float:
    R_fm = float(r0 * A ** (1 / 3))
    return R_fm / (HBAR_C_MEV_FM / ME_MEV)


def momentum(W: Real, rest_mass: float = 1.0) -> np.ndarray:
    """
    Helper function to get momentum value from energy.

    Args:
        W: Array of full particle energy values in natural units
        rest_mass: rest mass of the particle in natural units (m_e = 1 by default for electron)

    Returns:
        Array of particle momentum values.
    """
    return np.sqrt(W**2 - rest_mass**2)


def neutrino_energy(W: Real, W0: Real) -> Real:
    return W0 - W


def W_to_T(W: Real) -> np.ndarray:
    return np.asarray((W - 1.0) * ME_MEV, dtype=np.float64)


def T_to_W(T_MeV: Real) -> np.ndarray:
    return np.asarray(T_MeV / ME_MEV + 1.0, dtype=np.float64)
