# components/radiative.py
import logging

import numpy as np
from scipy.special import spence

from beta_spectrum.base import SpectrumComponent
from beta_spectrum.constants import (
    MP_MEV,
    ME_MEV,
    ALPHA,
    HBAR_C_MEV_FM,
    EULER_GAMMA,
    R0_FM,
)


class RadiativeCorrection(SpectrumComponent):
    """
    Outer radiative correction R(W, W0) = 1 + delta_1 + delta_2.

    Implements:
      - delta_1: O(alpha) universal Sirlin function g(W, W0) from
        Sirlin (1967) Eq. (20b), with endpoint soft-photon resummation
        per Sirlin (1987) Eq. (53).
      - delta_2: O(Z*alpha^2) finite nuclear size correction from
        Sirlin (1987) Eq. (1-10), using the modified Gaussian model.

    This is the INCLUSIVE correction suitable for:
      - Magnetic spectrometers
      - Detectors where bremsstrahlung is not detected
      - Precision at 1% level or worse

    The O(Z^2*alpha^3) correction is neglected for Z < 30 (spec Section 5).
    Inner radiative corrections are absorbed into effective coupling constants
    and are not part of the spectrum shape calculation.
    """

    def __init__(
        self,
        W0: float,
        Z: int = 0,
        A: int | None = None,
        use_endpoint_resummation: bool = True,
        delta_cut: float = 1e-3,
        logger: logging.Logger | None = None,
    ):
        """
        Args:
            W0: Endpoint energy in natural units m_e=c=h_bar=1
            Z: Nuclear charge of daughter nucleus
            A: Mass number (optional, affects nuclear model correction)
            use_endpoint_resummation: Apply Eq.(53) to handle ln(W0 - W) divergence
            delta_cut: Only apply resummation when (W0 - W) < delta_cut (spec Section 3.4)
            logger: Optional logger for debug/info output
        """
        super().__init__(logger=logger)
        self.W0 = W0
        self.Z = Z
        self.A = A
        self.use_endpoint_resummation = use_endpoint_resummation
        self.delta_cut = delta_cut
        self.m_p = MP_MEV / ME_MEV  # ~1836.15 (proton mass in m_e units)
        self.gamma_E = EULER_GAMMA  # Euler-Mascheroni constant

        # Precompute constant terms
        self._outer_normalization = ALPHA / (2.0 * np.pi)  # ~0.00116
        self._subtraction_term = 3.0 * np.log(self.m_p / (2.0 * self.W0))

        # Precompute O(Z*alpha^2) energy-independent nuclear part if A is given
        self._delta_F = self._compute_nuclear_model_correction(Z, A)

        if self._logger:
            self._logger.debug(
                "RadiativeCorrection: W0=%.4f, Z=%d, A=%s, delta_cut=%.2e, resummation=%s",
                W0,
                Z,
                A,
                delta_cut,
                use_endpoint_resummation,
            )

    def __call__(self, W: np.ndarray) -> np.ndarray:
        """
        Calculate R(W, W0) = 1 + delta_1(W, W0) + delta_2(Z, W).

        Args:
            W: Electron total energy array (in natural m_e=c=h_bar=1 units)

        Returns:
            Radiative correction factor array
        """
        if self._logger:
            self._logger.debug(
                "RadiativeCorrection: evaluating at %d points, W range=[%.4f, %.4f], W0=%.4f",
                len(W),
                W.min(),
                W.max(),
                self.W0,
            )

        W = np.asarray(W, dtype=float)

        # Avoid numerical issues exactly at endpoint W0
        eps = 1e-9
        mask_finite = W < (self.W0 - eps)

        delta_R = np.zeros_like(W, dtype=float)

        if np.any(mask_finite):
            W_finite = W[mask_finite]

            # O(alpha) correction (Sirlin 1967)
            delta_1 = self._delta_1(W_finite)

            # O(Z*alpha^2) correction (Sirlin 1987)
            delta_2 = self._delta_2(W_finite)

            delta_R_finite = delta_1 + delta_2
            delta_R[mask_finite] = delta_R_finite

        result = np.asarray(1.0 + delta_R, dtype=np.float64)

        if self._logger:
            self._logger.debug(
                "RadiativeCorrection: output range=[%.6e, %.6e]",
                result.min(),
                result.max(),
            )

        return result

    def _delta_1(self, W: np.ndarray) -> np.ndarray:
        """
        O(alpha) universal Sirlin function g(W, W0) per Sirlin (1967) Eq. (20b).

        g(W, W0) = 3*ln(m_p/(2*W0)) - 3/4
                   + 4*(tanh_inv_beta/beta - 1)*[(W0-W)/(3W) - 3/2 + ln(2*(W0-W))]
                   + (4/beta)*L(2*beta/(1+beta))
                   + (tanh_inv_beta/beta)*[2*(1+beta^2) + (W0-W)^2/(6*W^2) - 4*tanh_inv_beta]

        where L(x) = -Li_2(x) is the Spence function (dilogarithm).

        The correction is delta_1 = (alpha/2pi) * (g(W, W0) - subtraction_term).
        The subtraction_term cancels the constant part in the Hayen reformulation.
        """
        beta = self._safe_beta(W)
        tanh_inv_beta_beta = self._safe_tanh_inv_beta_over_beta(beta)

        # Dilogarithm term: L(2*beta/(1+beta)) = -Li_2(2*beta/(1+beta))
        arg = 2.0 * beta / (1.0 + beta)
        arg = np.minimum(arg, 1.0 - 1e-12)  # Avoid numerical problems at arg=1
        dilog_term = (4.0 / beta) * spence(1.0 - arg)

        # arctanh factor: (tanh^{-1}(beta)/beta) - 1
        atanh_factor = tanh_inv_beta_beta - 1.0

        # Energy difference
        delta_W = self.W0 - W
        delta_W = np.where(delta_W > 0, delta_W, 1e-12)

        if self.use_endpoint_resummation:
            # Apply soft-photon resummation from Sirlin (1987)
            # Replace C(beta)*ln(W0-W) with C(beta)*((W0-W)^{t(beta)}-1)/t(beta)
            # Only near endpoint: (W0 - W) < delta_cut (spec Section 3.4)
            t_beta = (2.0 * ALPHA / np.pi) * atanh_factor
            small_t = np.abs(t_beta) < 1e-6

            resummed_log = np.empty_like(delta_W)
            resummed_log[:] = np.log(delta_W)

            near_endpoint = delta_W < self.delta_cut
            if np.any(near_endpoint):
                resummed_log[near_endpoint] = (
                    np.power(delta_W[near_endpoint], t_beta[near_endpoint]) - 1.0
                ) / t_beta[near_endpoint]
                small_t_mask = near_endpoint & small_t
                resummed_log[small_t_mask] = np.log(delta_W[small_t_mask])

            log_term = np.log(2.0) + resummed_log
        else:
            log_term = np.log(2.0 * delta_W)

        # Energy-dependent term
        term2 = 4.0 * atanh_factor * (delta_W / (3.0 * W) - 1.5 + log_term)

        # Beta-dependent term (no divergence)
        term3 = tanh_inv_beta_beta * (
            2.0 * (1.0 + beta**2)
            + (delta_W**2) / (6.0 * W**2)
            - 4.0 * tanh_inv_beta_beta
        )

        # Leading constant
        const_term = -0.75

        g_func = const_term + dilog_term + term2 + term3

        return np.asarray(
            self._outer_normalization * (g_func - self._subtraction_term),
            dtype=np.float64,
        )

    def _delta_2(self, W: np.ndarray) -> np.ndarray:
        """
        O(Z*alpha^2) finite nuclear size correction per Sirlin (1987).

        delta_2(Z, W) = Z * alpha^2 * (Delta_1_4(W) + Delta_F)

        where:
          - Delta_1_4(W) = ln(m_p) - (5/12)*ln(2*W) + 43/18
            (energy-dependent relativistic part)
          - Delta_F(Z, A) = nuclear-structure-dependent part from modified
            Gaussian model (Eq. 10 of Sirlin 1987)

        For Z=0, this correction is zero.
        """
        if self.Z == 0:
            return np.zeros_like(W)

        # Energy-dependent part: Delta_1^0(W) + Delta_4(W)
        # From Sirlin 1987 relativistic limit
        Delta_1_4 = np.log(self.m_p) - (5.0 / 12.0) * np.log(2.0 * W) + 43.0 / 18.0

        # Total O(Z*alpha^2) correction
        delta_2 = self.Z * ALPHA**2 * (Delta_1_4 + self._delta_F)

        return np.asarray(delta_2, dtype=np.float64)

    def _compute_nuclear_model_correction(self, Z: int, A: int | None) -> float:
        """
        Compute the nuclear-structure-dependent part Delta_F of the O(Z*alpha^2) correction.

        Uses the modified Gaussian model from Sirlin (1987) Eq. (10):

        Delta_F = ln(Lambda/M) - kappa_1(Z)

        where:
          - Lambda = sqrt(6)/a, with a = r_0 * A^{1/3} (rms charge radius)
          - M = proton mass
          - kappa_1(Z) = 0.5 * [gamma_E + ln(3/(2*k^2)) + 2*alpha/(2+3*alpha)]
          - k^2 = (3/2) * (2+5*alpha)/(2+3*alpha)
          - alpha = (Z-2)/3
          - r_0 = 1.2 fm (Wilkinson)

        If A is not specified, returns 0 (no nuclear model correction).
        """
        if A is None:
            return 0.0

        # Nuclear radius from Wilkinson's r_0
        a_fm = R0_FM * A ** (1.0 / 3.0)  # rms charge radius in fm
        a_me = a_fm / (HBAR_C_MEV_FM / ME_MEV)  # convert to m_e^{-1} units

        # Lambda = sqrt(6)/a
        Lambda = np.sqrt(6.0) / a_me

        # M/Lambda ratio for the log term
        M_over_Lambda = self.m_p / Lambda

        # kappa_1(Z) from Sirlin 1987
        alpha_param = (Z - 2.0) / 3.0
        k_sq = 1.5 * (2.0 + 5.0 * alpha_param) / (2.0 + 3.0 * alpha_param)
        kappa_1 = 0.5 * (
            self.gamma_E
            + np.log(3.0 / (2.0 * k_sq))
            + 2.0 * ALPHA / (2.0 + 3.0 * ALPHA)
        )

        delta_F = float(np.log(M_over_Lambda) - kappa_1)
        return delta_F

    def _safe_beta(self, W: np.ndarray) -> np.ndarray:
        """
        Compute beta = p/W = sqrt(1 - 1/W^2) with safe handling near threshold.

        At W=1 (threshold), beta=0. We use a floor to avoid division by zero.
        """
        beta = np.sqrt(np.maximum(1.0 - 1.0 / W**2, 0.0))
        beta = np.maximum(beta, 1e-12)
        return beta

    def _safe_tanh_inv_beta_over_beta(self, beta: np.ndarray) -> np.ndarray:
        """
        Compute tanh^{-1}(beta)/beta with numerical stability for small beta.

        For beta < 1e-6, use Taylor expansion:
          tanh^{-1}(beta)/beta = 1 + beta^2/3 + beta^4/5 + ...

        For beta >= 1e-6, use the direct formula:
          tanh^{-1}(beta)/beta = 0.5 * ln((1+beta)/(1-beta)) / beta

        Per spec Section 3.3, this avoids division by zero and loss of precision.
        """
        result = np.empty_like(beta)
        small_beta = beta < 1e-6

        # Taylor series for small beta: 1 + beta^2/3 + beta^4/5 + beta^6/7
        beta_sq = beta**2
        result[small_beta] = (
            1.0
            + beta_sq[small_beta] / 3.0
            + beta_sq[small_beta] ** 2 / 5.0
            + beta_sq[small_beta] ** 3 / 7.0
        )

        # Direct computation for larger beta
        not_small = ~small_beta
        if np.any(not_small):
            result[not_small] = (
                0.5
                * np.log((1.0 + beta[not_small]) / (1.0 - beta[not_small]))
                / beta[not_small]
            )

        return result
