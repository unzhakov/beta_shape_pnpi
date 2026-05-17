"""
Theoretical model for fitting — wraps BetaSpectrum + DetectorResponse.

Provides a callable interface that evaluates the theoretical spectrum,
optionally convolved with detector response, on an arbitrary energy grid.

Usage:
    model = SpectrumModel(beta_spectrum, detector_response, channel_energies)
    y_model = model(channel_energies, normalization, background)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from beta_spectrum.spectrum import BetaSpectrum
    from beta_spectrum.components.detector_response import DetectorResponse
    from exp_data.spectrum import ExpSpectrum


class SpectrumModel:
    """Theoretical beta spectrum model for fitting to experimental data.

    Wraps a :class:`beta_spectrum.BetaSpectrum` and optionally a
    :class:`beta_spectrum.DetectorResponse` to produce a callable model
    that can be fitted against experimental spectra.

    The model computes:
        y(E) = A · [S(E) ⊗ R(E)] + B

    where S(E) is the theoretical spectrum, R(E) is the detector response,
    ⊗ denotes convolution, A is a normalization factor, and B is a constant
    background level.

    Parameters
    ----------
    spectrum : BetaSpectrum
        Theoretical spectrum calculator.
    detector_response : DetectorResponse, optional
        Detector response function. If None, no convolution is applied.
    channel_energies : np.ndarray, optional
        Energy values for detector channels (keV). Used to set up
        the response if detector_response is None.
    """

    def __init__(
        self,
        spectrum: "BetaSpectrum",
        detector_response: Optional["DetectorResponse"] = None,
        channel_energies: Optional[np.ndarray] = None,
    ):
        self.spectrum = spectrum
        self.detector_response = detector_response
        self.channel_energies = channel_energies

        # Cache for pre-computed unconvolved spectrum
        self._cached_W: Optional[np.ndarray] = None
        self._cached_spectral: Optional[np.ndarray] = None

    def evaluate_unconvolved(
        self,
        energies_keV: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Evaluate the theoretical spectrum (unconvolved).

        Parameters
        ----------
        energies_keV : np.ndarray, optional
            Energy grid in keV. Uses the spectrum's internal grid if None.

        Returns
        -------
        np.ndarray
            Theoretical spectrum values.
        """
        if energies_keV is not None:
            W = self._keV_to_W(energies_keV)
            return self.spectrum(W)
        return self.spectrum(
            self.spectrum._W_grid
            if hasattr(self.spectrum, "_W_grid") and self.spectrum._W_grid is not None
            else np.array([])
        )

    def _keV_to_W(self, energies_keV: np.ndarray) -> np.ndarray:
        """Convert kinetic energy in keV to total energy in m_e units."""
        return energies_keV / 510.998950 + 1.0

    def convolve(self, energies_keV: np.ndarray) -> np.ndarray:
        """Convolve theoretical spectrum with detector response.

        Parameters
        ----------
        energies_keV : np.ndarray
            Energy grid in keV for the detector channels.

        Returns
        -------
        np.ndarray
            Convolved spectrum on the given energy grid.

        Raises
        ------
        RuntimeError
            If no detector response is available.
        """
        if self.detector_response is None:
            raise RuntimeError(
                "No detector response configured. "
                "Provide a DetectorResponse to SpectrumModel."
            )

        # The detector's channel_energies are in m_e units
        W = self.detector_response.channel_energies
        # Clip W to valid range (W >= 1.0) before spectrum evaluation
        W_safe = np.maximum(W, 1.0 + 1e-12)
        spectral = self.spectrum(W_safe)

        # Clip negative/NaN values from spectrum (numerical artifacts)
        spectral = np.maximum(
            np.nan_to_num(spectral, nan=0.0, posinf=0.0, neginf=0.0), 0.0
        )

        # Convolve: detector_response.convolve(W, spectrum, normalize)
        return self.detector_response.convolve(W, spectral, normalize=True)

    def evaluate(
        self,
        energies_keV: np.ndarray,
        normalization: float = 1.0,
        background: float = 0.0,
    ) -> np.ndarray:
        """Evaluate the full theoretical model (with optional convolution).

        Parameters
        ----------
        energies_keV : np.ndarray
            Energy grid in keV.
        normalization : float
            Multiplicative normalization factor.
        background : float
            Constant background level.

        Returns
        -------
        np.ndarray
            Model prediction: A · (S ⊗ R) + B
        """
        if self.detector_response is not None:
            spectral = self.convolve(energies_keV)
        else:
            spectral = self.spectrum(self._keV_to_W(energies_keV))

        return normalization * spectral + background

    def residual(
        self,
        params: np.ndarray,
        exp_energies: np.ndarray,
        exp_counts: np.ndarray,
        exp_errors: np.ndarray,
    ) -> np.ndarray:
        """Compute weighted residuals for fitting.

        Parameters
        ----------
        params : np.ndarray
            Fit parameters [normalization, background, ...].
        exp_energies : np.ndarray
            Experimental energy grid (keV).
        exp_counts : np.ndarray
            Experimental counts.
        exp_errors : np.ndarray
            Experimental uncertainties.

        Returns
        -------
        np.ndarray
            Weighted residuals: (exp - model) / σ
        """
        model_values = self.evaluate(exp_energies, params[0], params[1])
        return (exp_counts - model_values) / exp_errors

    @staticmethod
    def default_params(
        exp_spectrum: "ExpSpectrum",
    ) -> tuple[float, float]:
        """Estimate default fit parameters from experimental data.

        Parameters
        ----------
        exp_spectrum : ExpSpectrum
            Experimental spectrum.

        Returns
        -------
        (normalization, background)
            Initial parameter estimates.
        """
        # Background from low-energy sideband
        mask_low = exp_spectrum.energies < np.percentile(exp_spectrum.energies, 5)
        background = (
            float(np.mean(exp_spectrum.counts[mask_low])) if np.any(mask_low) else 0.0
        )

        # Normalization from total counts / integral of unconvolved model
        total_counts = float(np.sum(exp_spectrum.counts))
        unconvolved = exp_spectrum.counts - background
        norm = total_counts / max(
            float(np.trapezoid(unconvolved, exp_spectrum.energies)), 1e-10
        )

        return (norm, background)

    @staticmethod
    def default_bounds(
        n_channels: int,
    ) -> tuple[list[float], list[float]]:
        """Estimate default parameter bounds.

        Parameters
        ----------
        n_channels : int
            Number of energy channels.

        Returns
        -------
        (lower_bounds, upper_bounds)
            Parameter bounds for [normalization, background].
        """
        return (
            [0.01, 0.0],
            [100.0, 1e6],
        )
