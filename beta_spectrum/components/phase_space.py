# components/phase_space.py
import numpy as np
from beta_spectrum.base import SpectrumComponent
from beta_spectrum.utils import momentum


class PhaseSpace(SpectrumComponent):
    """
    3-body phase space calculation (assuming heavy nucleus and no recoil)
    """

    def __init__(
        self,
        W0: float,
        transition_type: str = "A",
        m_e: float = 1.0,
        m_nu: float = 0.0,
    ):
        """
        Create instance of phase space with given decay parameters.

        Args:

            W0: Endpoint energy in natural units m_e=c=h_bar=1
            transition_type: allowed or N-order forbidden, should be [A, F1, F1U, F2, F2U, F3, F3U, F4]
            m_e: 1st body mass (electron by default: m_e = 1)
            m_nu:2nd body mass (massless neutrino by default: m_nu = 0)
        """
        self.tr_type = transition_type
        if self.tr_type not in ["A", "F1", "F1U", "F2", "F2U", "F3", "F3U", "F4"]:
            raise ValueError(
                f"Cannot calculate phase space for forbiddance '{self.tr_type}'.\nShould be [A, F1, F1U, F2, F2U, F3, F3U, F4].\n"
            )
        self.W0 = W0
        self.m_e = m_e
        self.m_nu = m_nu

    def __call__(self, W: np.ndarray) -> np.ndarray:
        w_e = W
        p_e = momentum(W, self.m_e)
        w_nu = self.W0 - W
        p_nu = momentum(w_nu, self.m_nu)

        phase_space = p_e * w_e * p_nu * w_nu

        if self.tr_type in ["A", "F1"]:
            forbid_factor = 1

        elif self.tr_type in ["F1U", "F2"]:
            forbid_factor = p_nu**2 + p_e**2

        elif self.tr_type in ["F2U", "F3"]:
            forbid_factor = p_nu**4 + 3 / 10 * p_nu**2 * p_e**2 + p_e**4

        elif self.tr_type in ["F3U", "F4"]:
            forbid_factor = (
                p_nu**6 + 7 * p_nu**4 * p_e**2 + 7 * p_nu**2 * p_e**4 + p_e**6
            )

        return np.asarray(phase_space * forbid_factor, dtype=np.float64)
