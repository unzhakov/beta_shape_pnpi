# units.py
import numpy as np

ME_KEV = 510.99895000  # keV electron mass


def W_to_T(W: float | np.ndarray) -> np.ndarray:
    return (W - 1.0) * ME_KEV


def T_to_W(T_keV: float | np.ndarray) -> np.ndarray:
    return T_keV / ME_KEV + 1.0
