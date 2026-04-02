# utils.py
from .constants import HBAR_C_MEV_FM, ME_MEV


def nuclear_radius(A: float) -> float:
    r0_fm = 1.2  # fm
    R_fm = r0_fm * A ** (1 / 3)

    return R_fm / (HBAR_C_MEV_FM / ME_MEV)
