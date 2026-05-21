"""Unit tests for the LivetimeDetermination class."""

from __future__ import annotations

import numpy as np
import pytest

from exp_data.spectrum import ExpSpectrum
from exp_data.livetime import LivetimeDetermination, LivetimeResult


class TestLivetimeResult:
    def test_basic(self):
        result = LivetimeResult(
            live_time_fraction=0.95,
            real_time_sec=3600.0,
            pulser_frequency_hz=100.0,
            observed_pulser_rate=95.0,
            pulser_peak_center_keV=10.0,
            pulser_peak_sigma_keV=0.5,
            pulser_peak_chi2=1.2,
        )
        assert result.live_time_fraction == 0.95
        assert result.metadata == {}

    def test_metadata_populated(self):
        result = LivetimeResult(
            live_time_fraction=0.95,
            real_time_sec=3600.0,
            pulser_frequency_hz=100.0,
            observed_pulser_rate=95.0,
            pulser_peak_center_keV=10.0,
            pulser_peak_sigma_keV=0.5,
            pulser_peak_chi2=1.2,
            metadata={"pulser_peak_channel": 100, "pulser_fit_success": True},
        )
        assert result.metadata["pulser_peak_channel"] == 100


class TestLivetimeDetermination:
    def test_from_file_requires_positive_frequency(self):
        """Test that zero frequency raises ValueError."""
        n_channels = 4096
        energies = np.arange(n_channels, dtype=np.float64)
        counts = np.ones(n_channels) * 10
        spectrum = ExpSpectrum(energies=energies, counts=counts)

        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=0.0,
            real_time_sec=3600.0,
        )
        with pytest.raises(ValueError, match="positive"):
            determiner.run()

    def test_live_time_clamped_to_1(self):
        """Test that live-time fraction > 1.0 is clamped to 1.0."""
        n_channels = 4096
        energies = np.arange(n_channels, dtype=np.float64)
        counts = np.exp(-energies / 500) * 10  # background
        # Add prominent pulser peak at channel 100
        pulser_channel = 100
        pulser_peak = 500 * np.exp(-0.5 * ((energies - pulser_channel) / 2.0) ** 2)
        counts += pulser_peak.astype(float)
        spectrum = ExpSpectrum(energies=energies, counts=counts)

        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=10.0,  # low frequency -> observed rate > frequency
            real_time_sec=0.01,  # very short time -> high observed rate
        )
        result = determiner.run()
        assert result.live_time_fraction == 1.0

    def test_live_time_clamped_to_0(self):
        """Test that live-time fraction < 0 is clamped to 0."""
        n_channels = 4096
        energies = np.arange(n_channels, dtype=np.float64)
        counts = np.exp(-energies / 500) * 10  # background only, no pulser peak
        spectrum = ExpSpectrum(energies=energies, counts=counts)

        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=100.0,
            real_time_sec=3600.0,
        )
        result = determiner.run()
        # With no pulser peak, the fitter may find a tiny spurious peak
        # but the observed rate will be negligible -> live_time ~ 0
        assert result.live_time_fraction >= 0.0
        assert result.live_time_fraction <= 0.01

    def test_pulser_peak_detection(self):
        """Test that LivetimeDetermination correctly identifies a pulser peak."""
        n_channels = 4096
        energies = np.arange(n_channels, dtype=np.float64)
        counts = np.exp(-energies / 500) * 10  # background
        # Add pulser peak at channel 100 (10 keV)
        pulser_channel = 100
        pulser_peak = 500 * np.exp(-0.5 * ((energies - pulser_channel) / 2.0) ** 2)
        counts += pulser_peak.astype(float)

        spectrum = ExpSpectrum(
            energies=energies,
            counts=counts,
            source="pulser",
        )

        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=100.0,
            real_time_sec=3600.0,
        )
        result = determiner.run()

        assert result.live_time_fraction > 0
        assert result.live_time_fraction <= 1.0
        assert result.pulser_peak_center_keV == pytest.approx(pulser_channel, abs=5.0)

    def test_result_contains_pulser_info(self):
        """Test that result contains pulser peak info."""
        n_channels = 4096
        energies = np.arange(n_channels, dtype=np.float64)
        counts = np.exp(-energies / 500) * 10
        pulser_channel = 100
        pulser_peak = 500 * np.exp(-0.5 * ((energies - pulser_channel) / 2.0) ** 2)
        counts += pulser_peak.astype(float)

        spectrum = ExpSpectrum(energies=energies, counts=counts, source="pulser")

        determiner = LivetimeDetermination(
            spectrum,
            pulser_frequency_hz=100.0,
            real_time_sec=3600.0,
        )
        result = determiner.run()

        assert result.pulser_peak_center_keV > 0
        assert result.pulser_peak_sigma_keV > 0
        assert result.pulser_peak_chi2 > 0
        assert result.pulser_frequency_hz == 100.0
        assert result.real_time_sec == 3600.0
