import numpy as np
from scipy.special import spence

from beta_spectrum.base import SpectrumComponent
from beta_spectrum.constants import MP_MEV, ME_MEV, ALPHA


class RadiativeCorrection(SpectrumComponent):
    """
    Outer radiative correction R(W, W0) = 1 + delta(W, W0).

    Implements Eq. (47)-(50) from Hayen et al. 2017, with endpoint
    divergence handeled by the soft photon resummation of Eq. (53).

    This is the INCLUSIVE correction suitable for:
    - Magnetic spectrometers
    - Detectors where bremmstrahlung is not detected
    - Precision at 1% level or worse
    """

    def __init__(self, W0: float, use_endpoint_resummation: bool = True):
        """
        Args:
            W0: Endpoint energy in natural units m_e=c=h_bar=1
            use_endpoint_resummation: Apply eq. (53) to handle ln(W0 - W) divergence
        """
        self.W0 = W0
        self.use_endpoint_resummation = use_endpoint_resummation
        self.m_p = MP_MEV / ME_MEV  # ~1836.15

        # Precompute constant terms of Eq. (48)
        self._outer_normalization = ALPHA / (2.0 * np.pi)  # ~0.00116
        self._subtraction_term = 3.0 * np.log(self.m_p / (2.0 * self.W0))

    def __call__(self, W: np.ndarray) -> np.ndarray:
        """
        Calculate R(W, W0) = 1 + delta_R(w, W0)

        Args:
            W: Electron total energy array (in natural m_e=c=h_bar=1 units)

        Returns:
            Radiative correction factor array
        """
        W = np.asarray(W)

        # Avoid numerical issues exactly at enpoint W0
        eps = 1e-9
        mask_finite = W < (self.W0 - eps)

        delta_R = np.zeros_like(W, dtype=float)

        if np.any(mask_finite):
            W_finite = W[mask_finite]

            if self.use_endpoint_resummation:
                g_func = self._g_function_resummed(W_finite)
            else:
                g_func = self._g_function_standard(W_finite)

            delta_R_finite = self._outer_normalization * (
                g_func - self._subtraction_term
            )
            delta_R[mask_finite] = delta_R_finite

        # Apply full correction including inner and L factor
        return 1.0 + delta_R

    def _g_function_standard(self, W: np.ndarray) -> np.ndarray:
        """
        Standard delta_R from Eq. (50), without endpoint resummation.
        Has logarithmic divergence as W -> W0
        """
        beta = np.sqrt(1.0 - 1.0 / W**2)
        beta = np.maximum(beta, 1e-12)  # Avoid division by 0

        # Dilogarithm term
        arg = 2.0 * beta / (1.0 + beta)
        arg = np.minimum(arg, 1.0 - 1e-12)  # Avoid numerical problems at arg=1
        dilog_term = (4.0 / beta) * spence(arg)

        # arctanh factor
        atanh_beta = np.arctanh(beta)
        atanh_factor = (atanh_beta / beta) - 1.0

        # Energy difference
        delta_W = self.W0 - W
        delta_W = np.maximum(delta_W, 1e-12)

        # Second line Eq. (50)
        term2 = 4.0 * atanh_factor * (delta_W / (3.0 * W) - 1.5 + np.log(2.0 * delta_W))

        # Third line Eq. (50)
        term3 = (atanh_beta / beta) * (
            2.0 * (1.0 + beta**2) + (delta_W**2) / (6.0 * W**2) - 4.0 * atanh_beta
        )

        # Leading constants from Eq. (50)
        const_term = 3.0 * np.log(self.m_p) - 0.75

        return const_term + dilog_term + term2 + term3

    def _g_function_resummed(self, W: np.ndarray) -> np.ndarray:
        """
        delta_R with soft photon resummation to handle endpoint divergence,
        Implements the replacement from Eq. (53):

        t(beta) * ln(W0 - W) -> (W0 - W)^t(beta) - 1
        """
        beta = np.sqrt(1.0 - 1.0 / W**2)
        beta = np.maximum(beta, 1e-12)

        # t(beta) from Eq. (54)
        atanh_beta = np.arctanh(beta)
        t_beta = (2.0 * ALPHA / np.pi) * (atanh_beta / beta - 1.0)

        # Dilogarithm term (same as _delta_R_standard)
        arg = 2.0 * beta / (1.0 + beta)
        arg = np.minimum(arg, 1.0 - 1e-12)
        dilog_term = (4.0 / beta) * spence(arg)

        # arctanh factor
        atanh_factor = (atanh_beta / beta) - 1.0

        # Energy difference
        delta_W = self.W0 - W
        delta_W = np.maximum(delta_W, 1e-12)

        # Apply resummation from Eq. (53)
        # Replace divergent ln(W0 - W) term with (W0 - W)^t(beta) - 1
        resummed_log_term = (delta_W**t_beta) - 1.0

        # Term2 with resummed log (instead of ln(2*(W0 - W)))
        term2_resummed = (
            4.0
            * atanh_factor
            * (delta_W / (3.0 * W) - 1.5 + 0.5 * resummed_log_term + np.log(2.0))
        )

        # Term3 (no divergence)
        term3 = (atanh_beta / beta) * (
            2.0 * (1.0 + beta**2) + (delta_W**2) / (6.0 * W**2) - 4.0 * atanh_beta
        )

        # Leading constants
        const_term = 3.0 * np.log(self.m_p) - 0.75

        return const_term + dilog_term + term2_resummed + term3
