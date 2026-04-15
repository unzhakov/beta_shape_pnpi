import numpy as np

from beta_spectrum.base import SpectrumComponent
from beta_spectrum.utils import momentum
from beta_spectrum.constants import ALPHA


class ScreeningCorrection(SpectrumComponent):
    """
    Atomic screning correction S(Z, W).

    Implements:
        - Ratio method (non-screened to screened F(Z,W)*p*W) using shifted energy W -> W_tilde
        - Smooth switching function to regularize low-energy behaviour

    S(W) = 1 + f(W) * (S_raw(W) - 1
    """

    def __init__(
        self,
        fermi_function,
        V0: float = None,
        Ws: float = 1.01,
        Delta: float = 0.01,
        eps: float = 1e-6,
        C=0.015,
    ):
        """
        Parameters
        ----------
        Z : nuclear charge
        fermi_function : callable
            Instance of FermiFunction class
        V0 : screening potential (in m_e units)
            If None, uses simple scaling law
        Ws : switching midpoint
        Delta : switching width
        eps : numerical floor for W_tilde
        """
        self.fermi = fermi_function
        self.Z = fermi_function.Z
        self.Ws = Ws
        self.Delta = Delta
        self.eps = eps
        self._C = C
        self.V0 = V0 if V0 is not None else self._estimate_V0()

    # --------------------
    # Core evaluation
    # --------------------
    def __call__(self, W: np.ndarray) -> np.ndarray:
        W = np.asarray(W)

        # Energy shift by screening potential V0
        W_tilde = self._screened_energy(W)

        # Momenta
        p = momentum(W)
        p_tilde = momentum(W_tilde)

        # Fermi functions
        F = self.fermi(W)
        F_tilde = self.fermi(W_tilde)

        # Raw ratio model (for W >> V0)
        S_raw = (F_tilde / F) * (p_tilde * W_tilde) / (p * W)

        # Smooth switching
        f = self._switching_function(W)

        return 1.0 + f * (S_raw - 1.0)

    # --------------------
    # Helpers
    # --------------------
    def _screened_energy(self, W):
        """
        Apply energy shift W -> W_tilde with safe numerical floor
        """
        W_tilde = W - self.V0
        return np.maximum(W_tilde, 1.0 + self.eps)

    def _switching_function(self, W):
        """
        Logistic switching function
        """
        return 1.0 / (1.0 + np.exp((self.Ws - W) / self.Delta))

    def _estimate_V0(self):
        """
        Simple scaling estimate for screening potential.

        V0 ~ alpha * Z^(4/3) * C

        C ~ 0.01-0.02 (tunable constant)
        """
        return ALPHA * self.Z ** (4.0 / 3.0) * self._C
