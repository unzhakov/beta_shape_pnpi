"""
Simple peak fitters for calibration source analysis.

Provides Gaussian and Gaussian+background fitting for identifying
known X-ray and gamma-ray lines used in energy calibration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.optimize import least_squares


@dataclass
class PeakFitResult:
    """Result of a peak fit."""

    center_keV: float
    """Peak center energy in keV."""

    sigma_keV: float
    """Peak width (σ) in keV."""

    amplitude: float
    """Peak amplitude (counts at peak)."""

    background: float
    """Constant background level."""

    chi2: float
    """Reduced χ²."""

    n_dof: int
    """Number of degrees of freedom."""

    covariance: np.ndarray
    """Parameter covariance matrix.

    Parameter order: [amplitude, center, sigma, background]
    """

    success: bool
    """Whether the fit converged."""

    def resolution(self) -> float:
        """FWHM in keV."""
        return 2.35482 * self.sigma_keV

    def resolution_fraction(self) -> float:
        """FWHM / center as a fraction."""
        if self.center_keV > 0:
            return self.resolution() / self.center_keV
        return 0.0


class GaussianFitter:
    """Fit a Gaussian peak on top of constant background.

    Model:  f(x) = A · exp(-(x - μ)² / (2σ²)) + B

    Parameters
    ----------
    energies : np.ndarray
        Energy values (keV), shape (N,).
    counts : np.ndarray
        Measured counts, shape (N,).
    uncertainties : np.ndarray, optional
        Measurement uncertainties. If None, assumed to be 1.0.
    """

    def __init__(
        self,
        energies: np.ndarray,
        counts: np.ndarray,
        uncertainties: Optional[np.ndarray] = None,
    ):
        self.energies = np.asarray(energies, dtype=np.float64)
        self.counts = np.asarray(counts, dtype=np.float64)

        if uncertainties is not None:
            self.uncertainties = np.asarray(uncertainties, dtype=np.float64)
        else:
            self.uncertainties = np.ones_like(self.counts)

        # Filter out zero-uncertainty points
        valid = self.uncertainties > 0
        self.energies = self.energies[valid]
        self.counts = self.counts[valid]
        self.uncertainties = self.uncertainties[valid]

    def _gaussian_background(
        self, params: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        """Gaussian + constant background model."""
        A, mu, sigma, B = params
        return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + B

    def _residuals(self, params: np.ndarray) -> np.ndarray:
        """Weighted residuals."""
        model = self._gaussian_background(params, self.energies)
        return (self.counts - model) / self.uncertainties

    def _chi2(self, params: np.ndarray) -> float:
        """Compute χ²."""
        residuals = self._residuals(params)
        return float(np.sum(residuals**2))

    def _compute_covariance(self, params: np.ndarray, jac: np.ndarray) -> np.ndarray:
        """Estimate covariance from Jacobian."""
        jtj = jac.T @ jac
        if np.linalg.cond(jtj) < 1e12:
            cov = np.linalg.inv(jtj)
            chi2 = self._chi2(params)
            ndof = max(len(self.energies) - 4, 1)
            cov = cov * chi2 / ndof
            return cov
        return np.diag([1e10, 1e6, 1e4, 1e6])

    def fit(
        self,
        guess: Optional[Tuple[float, float, float, float]] = None,
        bounds: Optional[Tuple[tuple, tuple]] = None,
    ) -> PeakFitResult:
        """Fit a Gaussian peak to the data.

        Parameters
        ----------
        guess : (amplitude, center, sigma, background), optional
            Initial parameter guess. Auto-estimated if None.
        bounds : ((A_min, μ_min, σ_min, B_min), (A_max, μ_max, σ_max, B_max)), optional
            Parameter bounds.

        Returns
        -------
        PeakFitResult
        """
        if guess is None:
            guess = self._estimate_guess()

        if bounds is None:
            bounds = (
                (0.0, min(self.energies), 0.1, 0.0),
                (np.inf, max(self.energies), 50.0, np.inf),
            )

        result = least_squares(
            self._residuals,
            x0=guess,
            bounds=bounds,
            method="trf",
            max_nfev=1000,
        )

        A, mu, sigma, B = result.x
        n_dof = max(len(self.energies) - 4, 1)
        chi2 = self._chi2(result.x)

        # Compute covariance from Jacobian
        _, jac = result.fun, result.jac
        cov = self._compute_covariance(result.x, jac) if jac is not None else np.eye(4)

        return PeakFitResult(
            center_keV=mu,
            sigma_keV=sigma,
            amplitude=A,
            background=B,
            chi2=chi2 / n_dof if n_dof > 0 else chi2,
            n_dof=n_dof,
            covariance=cov,
            success=result.success,
        )

    def _estimate_guess(self) -> Tuple[float, float, float, float]:
        """Auto-estimate initial parameters from data."""
        # Find the peak
        peak_idx = np.argmax(self.counts)
        center = self.energies[peak_idx]
        amplitude = self.counts[peak_idx]

        # Estimate sigma from FWHM
        half_max = amplitude / 2.0
        left_idx = np.argmax(self.counts[:peak_idx] < half_max)
        right_idx = np.argmin(self.counts[peak_idx:] < half_max) + peak_idx

        if left_idx < peak_idx and right_idx > peak_idx:
            fwhm = self.energies[right_idx] - self.energies[left_idx]
            sigma = fwhm / 2.35482
        else:
            sigma = 1.0  # Default width

        # Background estimate from average of edges
        n_edge = max(10, len(self.counts) // 20)
        background = np.mean(
            np.concatenate([self.counts[:n_edge], self.counts[-n_edge:]])
        )

        return (amplitude - background, center, max(sigma, 0.1), background)


class PeakFitter:
    """Find and fit multiple peaks in a spectrum.

    Parameters
    ----------
    energies : np.ndarray
        Energy values (keV).
    counts : np.ndarray
        Measured counts.
    uncertainties : np.ndarray, optional
        Measurement uncertainties.
    """

    def __init__(
        self,
        energies: np.ndarray,
        counts: np.ndarray,
        uncertainties: Optional[np.ndarray] = None,
    ):
        self.energies = np.asarray(energies, dtype=np.float64)
        self.counts = np.asarray(counts, dtype=np.float64)

        if uncertainties is not None:
            self.uncertainties = np.asarray(uncertainties, dtype=np.float64)
        else:
            self.uncertainties = np.sqrt(np.maximum(self.counts, 1.0))

    def find_peaks(
        self,
        min_height: float = 10.0,
        min_distance_keV: float = 2.0,
        smoothing_window: int = 5,
    ) -> np.ndarray:
        """Find candidate peak positions.

        Parameters
        ----------
        min_height : float
            Minimum peak height above background (counts).
        min_distance_keV : float
            Minimum separation between peaks in keV.
        smoothing_window : int
            Smoothing window for peak detection (odd number).

        Returns
        -------
        np.ndarray
            Energy values of detected peak centers.
        """
        from scipy.signal import find_peaks

        # Estimate background by smoothing
        window = smoothing_window if smoothing_window % 2 == 1 else smoothing_window + 1
        if len(self.counts) > window:
            from scipy.ndimage import uniform_filter1d

            background_est = uniform_filter1d(self.counts, size=window)
        else:
            background_est = np.full_like(self.counts, np.mean(self.counts))

        heights = self.counts - background_est
        indices, _ = find_peaks(
            heights, height=min_height, distance=int(min_distance_keV / max(1.0, np.diff(self.energies).min()))
        )

        return self.energies[indices]

    def fit_peaks(
        self,
        peak_energies: Optional[np.ndarray] = None,
        known_energies: Optional[np.ndarray] = None,
    ) -> list[PeakFitResult]:
        """Fit peaks at given or auto-detected positions.

        Parameters
        ----------
        peak_energies : np.ndarray, optional
            Explicit peak center energies to fit.
        known_energies : np.ndarray, optional
            Known calibration line energies — auto-detect peaks near these.

        Returns
        -------
        list[PeakFitResult]
        """
        if known_energies is not None:
            # Find peaks near known energies
            peak_energies = np.array([
                self._nearest_channel(e) for e in known_energies
            ])

        if peak_energies is None:
            peak_energies = self.find_peaks()

        results = []
        for center in peak_energies:
            fitter = GaussianFitter(self.energies, self.counts, self.uncertainties)
            result = fitter.fit(guess=(10.0, center, 1.0, 0.0))
            if result.success:
                results.append(result)

        return results

    def _nearest_channel(self, target_keV: float) -> float:
        """Find the energy channel closest to target_keV."""
        idx = np.argmin(np.abs(self.energies - target_keV))
        return self.energies[idx]
