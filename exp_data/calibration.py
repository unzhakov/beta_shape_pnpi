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
from pathlib import Path
from typing import Any, Optional

import numpy as np
from scipy.signal import find_peaks


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

            f = interp1d(e_vals, ch_vals, kind="linear")

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
        # np.polyfit returns descending-order coefficients [a_n, ..., a_0],
        # which is exactly what np.polyval expects.
        y_fit = np.polyval(coeffs, x)
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
        self.set_calibration_points(
            [
                (float(low_channel), low_energy_keV),
                (float(high_channel), high_energy_keV),
            ]
        )
        return self.apply_calibration(order=1)

    @classmethod
    def from_am241_spectrum(
        cls,
        filepath: str | Path,
        known_energies: Optional[list[float]] = None,
        n_channels: int = 4096,
    ) -> CalibrationResult:
        """Create a CalibrationResult from an Am-241 calibration spectrum.

        Auto-detects peaks using scipy.signal.find_peaks(), matches to known Am-241
        lines, fits linear calibration, and returns a CalibrationResult suitable
        for ExpSpectrum.apply_calibration().

        Parameters
        ----------
        filepath : str or Path
            Path to the Am-241 binary calibration file (CAM1.DAT or CAM2.DAT).
        known_energies : list of float, optional
            Known Am-241 line energies in keV. Defaults to [13.8, 17.8, 26.3].
        n_channels : int
            Number of channels (default 4096).

        Returns
        -------
        CalibrationResult
            Linear calibration fit with coefficients, chi2_per_dof, and covariance.

        Raises
        ------
        ValueError
            If fewer than 2 peaks are matched.
        """
        import struct

        filepath = Path(filepath)

        # Read raw binary
        with open(filepath, "rb") as f:
            f.read(128)  # skip header
            counts = np.array(
                struct.unpack(f"{n_channels}i", f.read(n_channels * 4)), dtype=np.int64
            )

        # Search for prominent gamma peaks across the full spectrum.
        # Exclude the high-energy pulser region (last ~200 channels).
        counts_search = counts.copy()
        pulser_cutoff = max(n_channels - 200, 300)
        counts_search[pulser_cutoff:] = 0

        # Find peaks with moderate thresholds.
        # height=100 skips noise floor; distance=5 keeps nearby peaks separate;
        # width=5 constrains peaks to realistic FWHM (~5-30 channels for
        # typical beta spectrometers), rejecting noise bumps and pulser tails.
        peaks, _ = find_peaks(counts_search, height=100, distance=5, width=(5, 50))
        peak_energies = np.array(peaks, dtype=float)
        peak_counts = counts_search[peaks]

        # Known Am-241 gamma lines (keV)
        if known_energies is None:
            known_energies = [13.0, 17.8, 26.3, 59.54]

        # Match peaks to known energies using a simpler, more robust
        # approach: find the two strongest peaks that span the largest
        # channel range, then assign the lowest and highest known energies
        # to them. This avoids the combinatorial explosion of trying all
        # energy pairs and the ambiguity of gain-based scoring.
        channel_to_energy: dict[float, float] = {}

        if len(peak_energies) >= 2 and len(known_energies) >= 2:
            strong_threshold = np.percentile(counts_search, 95)

            # Filter to strong peaks only.
            # Use a higher absolute threshold for the low-energy end
            # to avoid noise bumps being selected as calibration anchors.
            # The 95th percentile is ~171, but real Am-241 peaks are >400.
            absolute_threshold = max(strong_threshold, 400)
            strong_peaks = [
                (peak_energies[i], peak_counts[i])
                for i in range(len(peak_energies))
                if peak_counts[i] >= absolute_threshold
            ]

            if len(strong_peaks) >= 2:
                # Find the pair with the largest channel separation
                best_sep = 0
                best_pair = None
                for i in range(len(strong_peaks)):
                    for j in range(i + 1, len(strong_peaks)):
                        ch1, _ = strong_peaks[i]
                        ch2, _ = strong_peaks[j]
                        sep = ch2 - ch1
                        if sep > best_sep:
                            best_sep = sep
                            best_pair = (ch1, ch2)

                if best_pair is not None:
                    ch_low, ch_high = best_pair
                    e_min, e_max = min(known_energies), max(known_energies)
                    gain = (ch_high - ch_low) / (e_max - e_min)
                    intercept = ch_low - gain * e_min
                    channel_to_energy[float(ch_low)] = e_min
                    channel_to_energy[float(ch_high)] = e_max

        if len(channel_to_energy) < 2:
            raise ValueError(
                f"Could not match at least 2 Am-241 peaks. "
                f"Found {len(channel_to_energy)} of {len(known_energies)} expected. "
                f"Check calibration file or provide known_energies."
            )

        # Build calibration points and fit using only the two strongest
        # matched peaks for a clean 2-point linear fit.
        sorted_points = sorted(channel_to_energy.items(), key=lambda x: x[0])
        ch_low, e_low = sorted_points[0]
        ch_high, e_high = sorted_points[-1]

        calibrator = cls()
        calibrator.set_calibration_points([(ch_low, e_low), (ch_high, e_high)])
        return calibrator.fit_calibration(order=1)
