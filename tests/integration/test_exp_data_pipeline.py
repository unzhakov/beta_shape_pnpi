"""Integration tests for the experimental data processing pipeline."""

from __future__ import annotations

import numpy as np
import pytest

from exp_data.spectrum import ExpSpectrum
from exp_data.calibration import CalibrationResult
from exp_data.pipeline import SpectrumPipeline
from exp_data.livetime import LivetimeDetermination


def _has_root_io():
    """Check if any ROOT I/O library is available."""
    try:
        import ROOT  # noqa: F401

        return True
    except ImportError:
        pass
    try:
        import uproot  # noqa: F401

        return True
    except ImportError:
        return False


# --- ExpSpectrum.apply_calibration tests ---


class TestApplyCalibration:
    def test_returns_new_spectrum(self, synthetic_spectrum, calibration_result):
        """Test that apply_calibration returns a new ExpSpectrum."""
        calibrated = synthetic_spectrum.apply_calibration(calibration_result)
        assert calibrated is not synthetic_spectrum
        assert calibrated.n_channels == synthetic_spectrum.n_channels

    def test_energies_recalibrated(self, synthetic_spectrum, calibration_result):
        """Test that energies are recalibrated via the calibration polynomial."""
        calibrated = synthetic_spectrum.apply_calibration(calibration_result)
        expected = calibration_result(synthetic_spectrum.energies)
        np.testing.assert_allclose(calibrated.energies, expected)

    def test_counts_preserved(self, synthetic_spectrum, calibration_result):
        """Test that counts are preserved (copied, not modified)."""
        calibrated = synthetic_spectrum.apply_calibration(calibration_result)
        np.testing.assert_array_equal(calibrated.counts, synthetic_spectrum.counts)

    def test_errors_preserved(self, synthetic_spectrum, calibration_result):
        """Test that errors are preserved."""
        calibrated = synthetic_spectrum.apply_calibration(calibration_result)
        np.testing.assert_array_equal(calibrated.errors, synthetic_spectrum.errors)

    def test_metadata_updated_with_calibration(
        self, synthetic_spectrum, calibration_result
    ):
        """Test that metadata includes calibration info."""
        calibrated = synthetic_spectrum.apply_calibration(calibration_result)
        assert "calibration" in calibrated.metadata
        assert calibrated.metadata["calibration"]["type"] == "linear"
        assert (
            calibrated.metadata["calibration"]["chi2_per_dof"]
            == calibration_result.chi2_per_dof
        )

    def test_scalar_fields_preserved(self, synthetic_spectrum, calibration_result):
        """Test that all scalar fields are preserved."""
        spectrum = ExpSpectrum(
            energies=synthetic_spectrum.energies,
            counts=synthetic_spectrum.counts,
            dead_time=0.05,
            live_time=3420.0,
            source="Tc99",
            run_id="run_042",
            date="2025-01-15",
        )
        calibrated = spectrum.apply_calibration(calibration_result)
        assert calibrated.dead_time == 0.05
        assert calibrated.live_time == 3420.0
        assert calibrated.source == "Tc99"
        assert calibrated.run_id == "run_042"
        assert calibrated.date == "2025-01-15"

    def test_metadata_merged_not_replaced(self, synthetic_spectrum, calibration_result):
        """Test that existing metadata is merged, not replaced."""
        spectrum = ExpSpectrum(
            energies=synthetic_spectrum.energies,
            counts=synthetic_spectrum.counts,
            metadata={"existing_key": "existing_value"},
        )
        calibrated = spectrum.apply_calibration(calibration_result)
        assert "existing_key" in calibrated.metadata
        assert calibrated.metadata["existing_key"] == "existing_value"
        assert "calibration" in calibrated.metadata


# --- LivetimeDetermination tests ---


class TestLivetimeDetermination:
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

    def test_invalid_pulser_frequency(self):
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


# --- Pipeline validation tests ---


class TestPipelineValidation:
    def test_calib_chi2_error(self):
        """Test that poor calibration fit (chi2 > 3) produces an error message."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.arange(100, dtype=np.float64),
            counts=np.ones(100) * 10,
        )
        bad_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=5.0,  # poor fit
            n_points=2,
            covariance=np.eye(2),
        )
        messages = pipeline._validate(spectrum, bad_calib, None)
        assert any("ERROR" in m for m in messages)

    def test_calib_chi2_warning(self):
        """Test that marginal calibration fit (1.5 < chi2 < 3) produces a warning."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.arange(100, dtype=np.float64),
            counts=np.ones(100) * 10,
        )
        marginal_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=2.0,  # marginal fit
            n_points=2,
            covariance=np.eye(2),
        )
        messages = pipeline._validate(spectrum, marginal_calib, None)
        assert any("WARNING" in m for m in messages)
        assert not any("ERROR" in m for m in messages)

    def test_insufficient_counts(self):
        """Test that low counts produce an error."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.arange(10, dtype=np.float64),
            counts=np.array([1.0] * 10),  # only 10 counts total
        )
        good_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=0.5,
            n_points=2,
            covariance=np.eye(2),
        )
        messages = pipeline._validate(spectrum, good_calib, None)
        assert any("ERROR" in m and "counts" in m for m in messages)

    def test_negative_energy(self):
        """Test that negative energies are caught after calibration."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        # Create a spectrum with valid positive energies first
        spectrum = ExpSpectrum(
            energies=np.array([1.0, 2.0, 3.0, 4.0]),
            counts=np.array([10.0] * 4),
        )
        # Manually set energies to negative (simulating bad calibration)
        spectrum.energies = np.array([-1.0, 0.0, 1.0, 2.0])
        good_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=0.5,
            n_points=2,
            covariance=np.eye(2),
        )
        messages = pipeline._validate(spectrum, good_calib, None)
        assert any("ERROR" in m and "negative" in m.lower() for m in messages)

    def test_non_monotonic_energies(self):
        """Test that non-monotonic energies are caught."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 2.0, 1.0, 3.0]),
            counts=np.array([10.0] * 4),
        )
        good_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=0.5,
            n_points=2,
            covariance=np.eye(2),
        )
        messages = pipeline._validate(spectrum, good_calib, None)
        assert any("ERROR" in m and "monoton" in m for m in messages)

    def test_low_live_time_fraction_warning(self):
        """Test that low live-time fraction produces a warning."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.arange(100, dtype=np.float64),
            counts=np.ones(100) * 10,
        )
        good_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=0.5,
            n_points=2,
            covariance=np.eye(2),
        )
        bad_livetime = type(
            "MockLivetime",
            (),
            {
                "live_time_fraction": 0.3,
            },
        )()
        messages = pipeline._validate(spectrum, good_calib, bad_livetime)
        assert any("WARNING" in m and "live-time" in m.lower() for m in messages)

    def test_all_checks_pass(self):
        """Test that good data produces no errors."""
        pipeline = SpectrumPipeline(calibration_path="/dev/null")
        spectrum = ExpSpectrum(
            energies=np.arange(100, dtype=np.float64),
            counts=np.ones(100) * 100,
        )
        good_calib = CalibrationResult(
            coefficients=np.array([0.1, 0.4]),
            order=1,
            chi2_per_dof=0.5,
            n_points=2,
            covariance=np.eye(2),
        )
        good_livetime = type(
            "MockLivetime",
            (),
            {
                "live_time_fraction": 0.95,
            },
        )()
        messages = pipeline._validate(spectrum, good_calib, good_livetime)
        assert not any("ERROR" in m for m in messages)


# --- ROOT I/O integration tests ---


class TestRootIOIntegration:
    @pytest.mark.skipif(
        not _has_root_io(),
        reason="Neither PyROOT nor uproot available",
    )
    def test_round_trip(self, full_spectrum, tmp_path):
        """Test write and read round-trip preserves all fields."""
        from exp_data.root_io import write_spectrum, read_spectrum

        root_path = tmp_path / "test.root"
        write_spectrum(str(root_path), full_spectrum)

        read_back = read_spectrum(str(root_path))

        np.testing.assert_array_equal(read_back.energies, full_spectrum.energies)
        np.testing.assert_array_equal(read_back.counts, full_spectrum.counts)
        assert read_back.dead_time == full_spectrum.dead_time
        assert read_back.source == full_spectrum.source
        assert read_back.run_id == full_spectrum.run_id
        assert read_back.date == full_spectrum.date

    @pytest.mark.skipif(
        not _has_root_io(),
        reason="Neither PyROOT nor uproot available",
    )
    def test_metadata_round_trip(self, full_spectrum, tmp_path):
        """Test that metadata survives round-trip."""
        from exp_data.root_io import write_spectrum, read_spectrum

        root_path = tmp_path / "test.root"
        write_spectrum(str(root_path), full_spectrum)
        read_back = read_spectrum(str(root_path))

        assert "calibration" in read_back.metadata
        assert read_back.metadata["calibration"]["type"] == "linear"
        assert "livetime" in read_back.metadata
        assert read_back.metadata["livetime"]["fraction"] == 0.95


def _has_root_io():
    """Check if any ROOT I/O library is available."""
    try:
        import ROOT  # noqa: F401

        return True
    except ImportError:
        pass
    try:
        import uproot  # noqa: F401

        return True
    except ImportError:
        return False
