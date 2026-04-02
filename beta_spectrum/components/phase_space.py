# components/phase_space.py
import numpy as np
from ..base import SpectrumComponent


class PhaseSpace(SpectrumComponent):
    def __init__(self, W0: float, m_nu: float = 0.0):
        self.W0 = W0
        self.m_nu = m_nu

    def __call__(self, W: np.ndarray) -> np.ndarray:
        p = np.sqrt(W**2 - 1.0)
        Wv = self.W0 - W

        # neutrino mass generalization
        if self.m_nu == 0:
            return p * W * Wv**2
        else:
            return p * W * Wv * np.sqrt(Wv**2 - self.m_nu**2)
