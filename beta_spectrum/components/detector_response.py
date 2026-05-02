# components/detector_response.py
"""
Detector response modeling for beta spectrometry.

This module provides tools to convolve theoretical beta spectra with
realistic detector response functions, bridging the gap between
calculated spectra and measured data.

Two modes of operation:

1. **Tabulated response** — receive pre-computed response matrices
   (e.g. from GEANT4, EGSNRC, or dedicated calibration measurements).
   The response matrix has shape (N_calib_energies, N_channels), where
   each row is the detector response to a mono-energetic electron beam
   at the corresponding calibration energy.

2. **Analytical model** — parameterized Gaussian core + exponential tail,
   suitable for quick estimates and when no tabulated response is available.

The convolution operation:

    M_j = Σ_i Φ(W_i) · R(E_j; W_i) · ΔW_i

where Φ(W_i) is the theoretical differential spectrum at energy W_i,
R(E_j; W_i) is the probability that an electron of true energy W_i
deposits energy E_j in the detector channel j, and ΔW_i is the bin width.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d  # type: ignore[import-untyped]
from scipy.special import erfc  # type: ignore[import-untyped]
from typing import Optional

from beta_spectrum.constants import ME_MEV


class DetectorResponse:
    """
    Detector response function for beta spectrometry.

    Models the response of a semiconductor detector (e.g. PIPS, MMC)
    to mono-energetic electrons. Supports both tabulated responses
    (from MC simulations or calibration) and analytical parameterized models.

    Parameters
    ----------
    channel_energies : np.ndarray
        Energy values for each detector channel, shape (N_channels,).
        Typically a linear or log-spaced grid from 0 to Q-value.
    response_matrix : np.ndarray, optional
        Tabulated response matrix, shape (N_calib_energies, N_channels).
        If provided, analytical model parameters are ignored.
    calibration_energies : np.ndarray, optional
        Calibration energies corresponding to rows of response_matrix.
        Required if response_matrix is provided.
    model : str
        Analytical model to use when no tabulated response is given.
        Options: "gaussian", "gaussian_tail", "tikhonov".
    sigma_a : float
        Energy-independent contribution to resolution σ (in m_e units).
    sigma_b : float
        Coefficient for √E contribution to resolution σ (in m_e^(-1/2) units).
    tail_fraction : float
        Fraction of events in the exponential tail (0.0–1.0).
    tau : float
        Decay constant of exponential tail (in m_e units). Smaller = steeper tail.
    fano_factor : float
        Fano factor for energy resolution (typical: 0.12 for Si, 0.0 for MMC).
    """

    def __init__(
        self,
        channel_energies: np.ndarray,
        response_matrix: Optional[np.ndarray] = None,
        calibration_energies: Optional[np.ndarray] = None,
        model: str = "gaussian_tail",
        sigma_a: float = 0.0,
        sigma_b: float = 0.0,
        tail_fraction: float = 0.0,
        tau: float = 0.01,
        fano_factor: float = 0.12,
    ):
        self.channel_energies = np.asarray(channel_energies, dtype=float)
        self.n_channels = len(self.channel_energies)

        if response_matrix is not None and calibration_energies is not None:
            self._mode = "tabulated"
            self.response_matrix = np.asarray(response_matrix, dtype=float)
            self.calibration_energies = np.asarray(calibration_energies, dtype=float)
            self._validate_tabulated()
            self._interp = self._build_interpolator()
        else:
            self._mode = "analytical"
            self.model = model
            self.sigma_a = sigma_a
            self.sigma_b = sigma_b
            self.tail_fraction = tail_fraction
            self.tau = tau
            self.fano_factor = fano_factor

    def _validate_tabulated(self) -> None:
        """Validate tabulated response data."""
        if self.response_matrix.ndim != 2:
            raise ValueError(
                f"response_matrix must be 2D, got {self.response_matrix.ndim}D"
            )
        if self.response_matrix.shape[1] != self.n_channels:
            raise ValueError(
                f"response_matrix columns ({self.response_matrix.shape[1]}) "
                f"must match channel_energies length ({self.n_channels})"
            )
        if len(self.calibration_energies) != self.response_matrix.shape[0]:
            raise ValueError(
                f"calibration_energies length ({len(self.calibration_energies)}) "
                f"must match response_matrix rows ({self.response_matrix.shape[0]})"
            )
        if not np.all(np.diff(self.calibration_energies) > 0):
            raise ValueError("calibration_energies must be strictly increasing")
        if not np.all(np.diff(self.channel_energies) > 0):
            raise ValueError("channel_energies must be strictly increasing")

    def _build_interpolator(self) -> RegularGridInterpolator:
        """Build interpolation function for tabulated response."""
        return RegularGridInterpolator(
            (self.calibration_energies, self.channel_energies),
            self.response_matrix,
            method="linear",
            bounds_error=False,
            fill_value=0.0,
        )

    def get_response(self, energy: float) -> np.ndarray:
        """
        Get detector response array for a mono-energetic electron at given energy.

        Parameters
        ----------
        energy : float
            True electron energy in m_e units.

        Returns
        -------
        np.ndarray
            Response array of length N_channels.
        """
        if self._mode == "tabulated":
            # Interpolate from tabulated response
            resp = np.zeros(self.n_channels, dtype=float)
            for j in range(self.n_channels):
                f = interp1d(
                    self.calibration_energies,
                    self.response_matrix[:, j],
                    kind="linear",
                    bounds_error=False,
                    fill_value=0.0,
                )
                resp[j] = float(f(energy))
            return resp
        return self._analytical_response(energy)

    def _analytical_response(self, energy: float) -> np.ndarray:
        """
        Compute analytical detector response for mono-energetic electron.

        Uses a parameterized model:
        - Gaussian core with energy-dependent resolution
        - Optional exponential tail on the low-energy side
        """
        sigma = self._resolution_sigma(energy)

        if self.model == "gaussian":
            return self._gaussian_response(energy, sigma)
        if self.model == "gaussian_tail":
            return self._gaussian_tail_response(energy, sigma)
        if self.model == "tikhonov":
            return self._tikhonov_response(energy, sigma)
        raise ValueError(f"Unknown model: {self.model}")

    def _resolution_sigma(self, energy: float) -> float:
        """
        Compute energy-dependent resolution σ.

        σ(E) = sqrt(σ₀² + (b·√E)²)
        where σ₀ = sigma_a (energy-independent term)
              b = sigma_b (√E-dependent term)

        For semiconductor detectors, can also include Fano factor contribution:
        σ_Fano = sqrt(F · ε · w) where F is Fano factor, ε is energy,
        w is mean ionization energy (~3.6 eV for Si).
        """
        E = max(energy, 1e-12)
        sigma_sq = self.sigma_a**2 + (self.sigma_b * np.sqrt(E)) ** 2
        if self.fano_factor > 0:
            # Fano contribution: σ = sqrt(F * E * w) in natural units
            w_me = 3.6e-3 / ME_MEV  # 3.6 eV in m_e units
            sigma_sq += self.fano_factor * E * w_me
        return float(np.sqrt(sigma_sq))

    def _gaussian_response(self, energy: float, sigma: float) -> np.ndarray:
        """Pure Gaussian response."""
        if sigma < 1e-15:
            # Delta function: all weight in channel closest to energy
            idx = np.argmin(np.abs(self.channel_energies - energy))
            resp = np.zeros(self.n_channels, dtype=float)
            resp[idx] = 1.0
            return resp

        delta_E = self.channel_energies - energy
        resp = np.exp(-0.5 * (delta_E / sigma) ** 2)
        resp /= sigma * np.sqrt(2.0 * np.pi)
        return np.asarray(resp, dtype=np.float64)

    def _gaussian_tail_response(self, energy: float, sigma: float) -> np.ndarray:
        """
        Gaussian core + exponential tail on the low-energy side.

        f(E) = (1 - f_tail) · G(E; E₀, σ)
               + f_tail · T(E; E₀, σ, τ)

        where T is the tail function:
        T(E) = (1/(2τ)) · exp(τ·(E - E₀)/2 + σ²·τ²/8) · erfc(
                   (E - E₀)/(√2·σ) + σ·τ/√2
               )
        for E < E₀, and zero for E > E₀.

        This is the standard model for Si detector response with charge
        collection inefficiencies causing low-energy tailing.
        """
        if sigma < 1e-15:
            idx = np.argmin(np.abs(self.channel_energies - energy))
            resp = np.zeros(self.n_channels, dtype=float)
            resp[idx] = 1.0
            return resp

        delta_E = self.channel_energies - energy
        core = np.exp(-0.5 * (delta_E / sigma) ** 2)
        core /= sigma * np.sqrt(2.0 * np.pi)

        if self.tail_fraction <= 0 or self.tau <= 0:
            return np.asarray(core, dtype=np.float64)

        # Exponential tail (Tikhonov function): defined for all E
        tail = np.zeros(self.n_channels, dtype=float)
        dE_all = delta_E
        arg = dE_all / (np.sqrt(2.0) * sigma) + sigma * self.tau / np.sqrt(2.0)
        tail = (
            1.0
            / (2.0 * self.tau)
            * np.exp(self.tau * dE_all / 2.0 + (sigma * self.tau) ** 2 / 8.0)
            * erfc(arg)
        )

        resp = (1.0 - self.tail_fraction) * core + self.tail_fraction * tail

        # Normalize to unit area
        dw = np.diff(self.channel_energies)
        dw = np.concatenate([dw, [dw[-1]]])
        resp *= dw
        area = float(np.sum(resp))
        resp /= dw
        if area > 0:
            resp /= area

        return np.asarray(resp, dtype=np.float64)

    def _tikhonov_response(self, energy: float, sigma: float) -> np.ndarray:
        """
        Full Tikhonov function for detector response.

        Combines Gaussian core with exponential tail in a single
        analytic expression. Standard model for semiconductor detectors.
        """
        return self._gaussian_tail_response(energy, sigma)

    def convolve(
        self,
        W: np.ndarray,
        spectrum_values: np.ndarray,
        normalize: bool = True,
    ) -> np.ndarray:
        """
        Convolve theoretical spectrum with detector response.

        Computes the predicted measured spectrum:
            M_j = Σ_i Φ(W_i) · R(E_j; W_i) · ΔW_i

        Parameters
        ----------
        W : np.ndarray
            Energy grid of theoretical spectrum, shape (N_bins,).
        spectrum_values : np.ndarray
            Differential spectrum dN/dW at energies W, shape (N_bins,).
        normalize : bool
            If True, normalize each response function to unit area.

        Returns
        -------
        np.ndarray
            Predicted measured spectrum, shape (N_channels,).
        """
        W = np.asarray(W, dtype=float)
        spectrum_values = np.asarray(spectrum_values, dtype=float)

        if W.ndim != 1 or spectrum_values.ndim != 1:
            raise ValueError("Inputs must be 1D arrays")
        if len(W) != len(spectrum_values):
            raise ValueError("W and spectrum_values must have same length")

        # Compute bin widths
        dw = self._compute_bin_widths(W)

        # Build response matrix: shape (N_bins, N_channels)
        response_matrix = np.zeros((len(W), self.n_channels))

        for i, w_i in enumerate(W):
            resp = self.get_response(w_i)
            response_matrix[i] = resp

        # Convolve: M = R^T · (Φ · ΔW)
        convolved = response_matrix.T @ (spectrum_values * dw)
        return np.asarray(convolved, dtype=np.float64)

    def convolve_batch(
        self,
        W: np.ndarray,
        spectrum_values: np.ndarray,
        normalize: bool = True,
    ) -> np.ndarray:
        """
        Convolve with optimized batch computation for analytical models.

        For analytical models, computes the full response matrix at once
        rather than looping over energies. Much faster for large spectra.
        """
        if self._mode != "analytical":
            return self.convolve(W, spectrum_values, normalize)

        W = np.asarray(W, dtype=float)
        spectrum_values = np.asarray(spectrum_values, dtype=float)

        dw = self._compute_bin_widths(W)

        # Build full response matrix at once: shape (N_bins, N_channels)
        response_matrix = np.zeros((len(W), self.n_channels))

        for i, w_i in enumerate(W):
            if spectrum_values[i] <= 0:
                continue
            resp = self._analytical_response(w_i)
            response_matrix[i] = resp

        return np.asarray(response_matrix.T @ (spectrum_values * dw), dtype=np.float64)

    @staticmethod
    def _compute_bin_widths(W: np.ndarray) -> np.ndarray:
        """Compute effective bin widths for arbitrary energy grid."""
        if len(W) == 1:
            return np.array([1.0])
        dw = np.diff(W)
        # Half-bin widths at each point
        widths = np.zeros(len(W))
        widths[0] = dw[0]
        widths[-1] = dw[-1]
        widths[1:-1] = 0.5 * (dw[:-1] + dw[1:])
        return widths

    @classmethod
    def from_gaussian_params(
        cls,
        channel_energy_range: tuple[float, float],
        n_channels: int = 4096,
        sigma_a: float = 0.0,
        sigma_b: float = 0.0,
        tail_fraction: float = 0.0,
        tau: float = 0.01,
        model: str = "gaussian_tail",
        fano_factor: float = 0.12,
    ) -> DetectorResponse:
        """
        Convenience factory: create DetectorResponse with analytical Gaussian model.

        Parameters
        ----------
        channel_energy_range : tuple
            (E_min, E_max) in m_e units.
        n_channels : int
            Number of detector channels (default 4096).
        sigma_a, sigma_b : float
            Resolution parameters: σ(E) = sqrt(σ_a² + (σ_b·√E)²).
        tail_fraction : float
            Fraction of events in exponential tail (0.0–1.0).
        tau : float
            Exponential tail decay constant.
        model : str
            "gaussian", "gaussian_tail", or "tikhonov".
        fano_factor : float
            Fano factor (0.12 for Si, ~0 for MMC).
        """
        channel_energies = np.linspace(
            channel_energy_range[0],
            channel_energy_range[1],
            n_channels,
        )
        return cls(
            channel_energies=channel_energies,
            model=model,
            sigma_a=sigma_a,
            sigma_b=sigma_b,
            tail_fraction=tail_fraction,
            tau=tau,
            fano_factor=fano_factor,
        )

    @classmethod
    def from_mc_simulation(
        cls,
        channel_energies: np.ndarray,
        response_matrix: np.ndarray,
        calibration_energies: np.ndarray,
    ) -> DetectorResponse:
        """
        Convenience factory: create DetectorResponse from MC simulation output.

        Parameters
        ----------
        channel_energies : np.ndarray
            Energy per channel, shape (N_channels,).
        response_matrix : np.ndarray
            Response matrix from MC, shape (N_calib_energies, N_channels).
            Each row i is the response to mono-energetic electrons at
            calibration_energies[i].
        calibration_energies : np.ndarray
            Calibration energies, shape (N_calib_energies,).
        """
        return cls(
            channel_energies=channel_energies,
            response_matrix=response_matrix,
            calibration_energies=calibration_energies,
        )

    def __repr__(self) -> str:
        mode_str = "tabulated" if self._mode == "tabulated" else "analytical"
        if self._mode == "tabulated":
            return (
                f"DetectorResponse(mode={mode_str}, "
                f"calib_energies={len(self.calibration_energies)}, "
                f"channels={self.n_channels})"
            )
        else:
            return (
                f"DetectorResponse(mode={mode_str}, model={self.model!r}, "
                f"sigma_a={self.sigma_a:.4f}, sigma_b={self.sigma_b:.4f}, "
                f"tail_frac={self.tail_fraction:.3f}, "
                f"channels={self.n_channels})"
            )
