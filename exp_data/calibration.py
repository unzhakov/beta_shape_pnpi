"""
Energy calibration for beta spectrometers.

Provides methods to calibrate detector channel numbers to physical
energy values using known spectral lines (X-rays, gamma-rays).

Supports:
- Linear calibration: E = a + b · channel
- Quadratic calibration: E = a + b · channel + c · channel²
- Multi-point calibration with peak fitting
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np



@dataclass
class CalibrationResult:
    """Result of an energy calibration fit."""

    coefficients: np.ndarray
    """Calibration coefficients [a, b, c, ...] depending on polynomial order."""

    order: int
    """Polynomial order (1 = linear, 2 = quadratic)."""

    chi2_per_dof: float
    """Reduced χ² of the calibration fit."""

    n_points: int
    """Number of calibration points used."""

    covariance: np.ndarray
    """Parameter covariance matrix."""

    def __call__(self, channel: float | np.ndarray) -> np.ndarray:
        """Convert channel number(s) to energy (keV)."""
        return np.polyval(self.coefficients, channel)

    def invert(self) -> np.poly1d:
        """Return polynomial that converts energy → channel."""
        if self.order == 1:
            # coefficients from polyfit: [slope, intercept]
            b, a = self.coefficients
            # E = a + b*ch => ch = (E - a) / b
            # np.poly1d([1/b, -a/b]) => (1/b)*E + (-a/b) = (E-a)/b ✓
            return np.poly1d([1.0 / b, -a / b])
        else:
            # For quadratic, compute inverse numerically via interpolation
            ch_vals = np.linspace(-10, 200, 1000)
            e_vals = np.polyval(self.coefficients[::-1], ch_vals)
            from scipy.interpolate import interp1d
            f = interp1d(e_vals, ch_vals, kind='linear')
            class _InversePoly:
                def __init__(self, fn: Any) -> None:
                    self._fn = fn
                def __call__(self, x: float | np.ndarray) -> np.ndarray:
                    return self._fn(x)
            return _InversePoly(f)

    def resolution_at(self, energy_keV: float, sigma_channel: float = 1.0) -> float:
        """Estimate energy resolution (σ) at a given energy.

        Parameters
        ----------
        energy_keV : float
            Energy in keV.
        sigma_channel : float
            Channel width in channel units (typically 1.0 for integer channels).

        Returns
        -------
        float
            Energy resolution σ in keV.
        """
        # dE/dch = derivative of calibration polynomial at energy_keV
        # coefficients are stored as [a_n, ..., a_1, a_0]
        deriv_coeffs = []
        coeffs = self.coefficients
        for i in range(len(coeffs) - 1):
            deriv_coeffs.append(coeffs[i] * (len(coeffs) - i - 1))
        if not deriv_coeffs:
            return abs(sigma_channel)
        dE_dch = np.polyval(deriv_coeffs, energy_keV)
        return abs(dE_dch) * sigma_channel


class EnergyCalibrator:
    """Energy calibration for detector spectra.

    Parameters
    ----------
    channel_energies : np.ndarray
        Energy values (keV) per channel, shape (N_channels,).
        If None, channels are assumed to be integer indices.
    """

    def __init__(self, channel_energies: Optional[np.ndarray] = None):
        self._channel_energies = channel_energies
        self._calibration: Optional[CalibrationResult] = None
        self._calibration_points: Optional[list[tuple[float, float]]] = None

    @property
    def is_calibrated(self) -> bool:
        """Whether a calibration has been applied."""
        return self._calibration is not None

    def add_calibration_point(
        self,
        channel: int | float,
        energy_keV: float,
    ) -> None:
        """Add a calibration point (known line energy → channel number).

        Parameters
        ----------
        channel : int or float
            Detector channel number.
        energy_keV : float
            Known energy of the calibration line (keV).
        """
        if self._calibration_points is None:
            self._calibration_points = []
        self._calibration_points.append((float(channel), float(energy_keV)))

    def set_calibration_points(self, points: list[tuple[float, float]]) -> None:
        """Set calibration points directly.

        Parameters
        ----------
        points : list of (channel, energy_keV)
            Known calibration line positions.
        """
        self._calibration_points = [(float(c), float(e)) for c, e in points]

    def fit_calibration(
        self,
        order: int = 1,
        channel_energies: Optional[np.ndarray] = None,
    ) -> CalibrationResult:
        """Fit calibration polynomial using known points.

        Parameters
        ----------
        order : int
            Polynomial order (1 = linear, 2 = quadratic).
        channel_energies : np.ndarray, optional
            Full channel-to-energy mapping. If provided, fits to
            all points; otherwise fits only to added calibration points.

        Returns
        -------
        CalibrationResult
        """
        if channel_energies is not None:
            # Fit to full mapping
            channels = np.arange(len(channel_energies), dtype=np.float64)
            x = channels
            y = np.asarray(channel_energies, dtype=np.float64)
        elif self._calibration_points is not None:
            x = np.array([p[0] for p in self._calibration_points])
            y = np.array([p[1] for p in self._calibration_points])
        else:
            raise ValueError(
                "No calibration data available. Use add_calibration_point() "
                "or pass channel_energies."
            )

        # Fit polynomial: E(ch) = a₀ + a₁·ch + a₂·ch² + ...
        try:
            coeffs, cov = np.polyfit(x, y, order, cov=True)
        except ValueError:
            # Too few points for covariance estimation
            coeffs = np.polyfit(x, y, order)
            cov = np.eye(order + 1) * 1e-10

        # Compute χ²
        y_fit = np.polyval(coeffs[::-1], x)
        # Assume uniform uncertainty of 0.1 keV for calibration points
        uncertainties = np.full_like(x, 0.1)
        chi2 = float(np.sum(((y - y_fit) / uncertainties) ** 2))
        n_dof = max(len(x) - order - 1, 1)

        return CalibrationResult(
            coefficients=coeffs,
            order=order,
            chi2_per_dof=chi2 / n_dof,
            n_points=len(x),
            covariance=cov,
        )

    def apply_calibration(
        self,
        order: int = 1,
        channel_energies: Optional[np.ndarray] = None,
    ) -> CalibrationResult:
        """Fit and store calibration for future use.

        See :meth:`fit_calibration` for parameters.
        """
        self._calibration = self.fit_calibration(order, channel_energies)
        return self._calibration

    def channel_to_energy(self, channel: float | np.ndarray) -> np.ndarray:
        """Convert channel number to energy using current calibration.

        Parameters
        ----------
        channel : float or np.ndarray
            Channel number(s).

        Returns
        -------
        np.ndarray
            Energy in keV.

        Raises
        ------
        RuntimeError
            If no calibration has been applied.
        """
        if self._calibration is None:
            raise RuntimeError(
                "No calibration applied. Call apply_calibration() first."
            )
        return self._calibration(channel)

    def energy_to_channel(self, energy_keV: float | np.ndarray) -> np.ndarray:
        """Convert energy to channel using current calibration.

        Parameters
        ----------
        energy_keV : float or np.ndarray
            Energy in keV.

        Returns
        -------
        np.ndarray
            Channel number(s).

        Raises
        ------
        RuntimeError
            If no calibration has been applied.
        """
        if self._calibration is None:
            raise RuntimeError(
                "No calibration applied. Call apply_calibration() first."
            )
        inv = self._calibration.invert()
        return inv(energy_keV)

    def get_linear_calibration(
        self,
        low_energy_keV: float,
        high_energy_keV: float,
        low_channel: int,
        high_channel: int,
    ) -> CalibrationResult:
        """Create a linear calibration from two known points.

        Parameters
        ----------
        low_energy_keV : float
            Energy of first calibration line (keV).
        high_energy_keV : float
            Energy of second calibration line (keV).
        low_channel : int
            Channel number of first line.
        high_channel : int
            Channel number of second line.

        Returns
        -------
        CalibrationResult
        """
        self.set_calibration_points([
            (float(low_channel), low_energy_keV),
            (float(high_channel), high_energy_keV),
        ])
        return self.apply_calibration(order=1)
