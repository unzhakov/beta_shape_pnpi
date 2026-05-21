"""
Live-time determination from pulser calibration spectra.

Uses Gaussian peak fitting to identify the pulser peak in a calibration
spectrum and computes the live-time fraction from the observed peak
count rate compared to the known injection frequency.

Usage:
    from exp_data.livetime import LivetimeDetermination

    # From an existing ExpSpectrum
    determiner = LivetimeDetermination(
        spectrum=pulser_spectrum,
        pulser_frequency_hz=100.0,
        real_time_sec=3600.0,
    )
    result = determiner.run()

    # Or directly from a binary file
    result = LivetimeDetermination.from_file(
        filepath="data/raw/CAM1.DAT",
        pulser_frequency_hz=100.0,
        real_time_sec=3600.0,
    )
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from exp_data.spectrum import ExpSpectrum
from exp_data.raw_reader import read_raw_spectrum
from exp_data.fitters import GaussianFitter


@dataclass
class LivetimeResult:
    """Result of a live-time determination."""

    live_time_fraction: float
    """Live-time fraction tau (0.0 to 1.0)."""

    real_time_sec: float
    """Total real acquisition time in seconds."""

    pulser_frequency_hz: float
    """Known pulser injection frequency in Hz."""

    observed_pulser_rate: float
    """Observed pulser peak count rate in Hz."""

    pulser_peak_center_keV: float
    """Energy of the pulser peak in keV (from Gaussian fit)."""

    pulser_peak_sigma_keV: float
    """Width of the pulser peak in keV."""

    pulser_peak_chi2: float
    """Reduced chi-squared of the pulser peak fit."""

    metadata: dict = field(default_factory=dict)
    """Additional metadata from the determination."""


class LivetimeDetermination:
    """Determine live-time fraction from a pulser calibration spectrum.

    The pulser injects known-frequency pulses into the MCA. The observed
    peak count rate in the pulser channel, compared to the known injection
    frequency, gives the live-time fraction:

        tau = observed_pulser_rate / expected_pulser_frequency

    Parameters
    ----------
    spectrum : ExpSpectrum
        Calibration spectrum with pulser input (energies in keV).
    pulser_frequency_hz : float
        Known pulser injection frequency in Hz.
    real_time_sec : float
        Total real acquisition time in seconds (from MCA metadata).
    pulser_channel_guess : float, optional
        Expected pulser peak channel number (for narrowing search).
        If None, searches the full spectrum.
    """

    def __init__(
        self,
        spectrum: ExpSpectrum,
        pulser_frequency_hz: float,
        real_time_sec: float,
        pulser_channel_guess: Optional[float] = None,
    ):
        self.spectrum = spectrum
        self.pulser_frequency_hz = pulser_frequency_hz
        self.real_time_sec = real_time_sec
        self.pulser_channel_guess = pulser_channel_guess

    def run(self) -> LivetimeResult:
        """Determine live-time fraction from the pulser peak.

        Returns
        -------
        LivetimeResult
            Live-time fraction and associated metadata.

        Raises
        ------
        ValueError
            If no pulser peak is found or fit fails.
        """
        energies = self.spectrum.energies
        counts = self.spectrum.counts
        uncertainties = (
            self.spectrum.errors
            if self.spectrum.errors is not None
            else np.sqrt(np.maximum(counts, 0))
        )

        # Narrow search around pulser channel guess if provided
        if self.pulser_channel_guess is not None:
            mask = np.abs(energies - self.pulser_channel_guess) < 5.0
            search_energies = energies[mask]
            search_counts = counts[mask]
            search_uncertainties = uncertainties[mask]
        else:
            search_energies = energies
            search_counts = counts
            search_uncertainties = uncertainties

        # Find the highest peak (pulser should be prominent)
        peak_idx = np.argmax(search_counts)
        peak_center = search_energies[peak_idx]

        # Fit Gaussian to pulser peak
        fitter = GaussianFitter(search_energies, search_counts, search_uncertainties)
        result = fitter.fit(
            guess=(search_counts[peak_idx], peak_center, 1.0, np.mean(search_counts))
        )

        if not result.success:
            raise ValueError(f"Pulser peak fit failed: {result}")

        # Compute observed rate from peak amplitude
        # Peak amplitude is in counts; divide by real_time to get rate
        observed_rate = (
            result.amplitude / self.real_time_sec if self.real_time_sec > 0 else 0.0
        )

        # Compute live-time fraction
        if self.pulser_frequency_hz <= 0:
            raise ValueError("Pulser frequency must be positive")
        live_time_fraction = observed_rate / self.pulser_frequency_hz

        # Clamp to valid range [0, 1]
        live_time_fraction = max(0.0, min(1.0, live_time_fraction))

        return LivetimeResult(
            live_time_fraction=live_time_fraction,
            real_time_sec=self.real_time_sec,
            pulser_frequency_hz=self.pulser_frequency_hz,
            observed_pulser_rate=observed_rate,
            pulser_peak_center_keV=result.center_keV,
            pulser_peak_sigma_keV=result.sigma_keV,
            pulser_peak_chi2=result.chi2,
            metadata={
                "pulser_peak_channel": int(peak_center),
                "pulser_fit_success": result.success,
            },
        )

    @staticmethod
    def from_file(
        filepath: str | Path,
        pulser_frequency_hz: float,
        real_time_sec: float,
        n_channels: int = 4096,
        run_id: Optional[str] = None,
        source: Optional[str] = None,
    ) -> LivetimeResult:
        """Convenience: read spectrum from file and determine live-time.

        Parameters
        ----------
        filepath : str or Path
            Path to the pulser calibration binary file.
        pulser_frequency_hz : float
            Known pulser injection frequency in Hz.
        real_time_sec : float
            Total real acquisition time in seconds.
        n_channels : int
            Number of channels.
        run_id : str, optional
            Run identifier.
        source : str, optional
            Source identifier (e.g., 'Am241', 'pulser').

        Returns
        -------
        LivetimeResult
        """
        spectrum = read_raw_spectrum(
            filepath,
            n_channels=n_channels,
            run_id=run_id,
            source=source,
            live_time_sec=real_time_sec,
            metadata={"operation": "livetime_calibration"},
        )
        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=pulser_frequency_hz,
            real_time_sec=real_time_sec,
        )
        return determiner.run()
