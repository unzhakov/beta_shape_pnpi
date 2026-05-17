"""Tests for exp_data.spectrum."""

import numpy as np
import pytest

from exp_data.spectrum import ExpSpectrum


class TestExpSpectrumBasic:
    """Basic ExpSpectrum construction and properties."""

    def test_construct_with_errors(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0, 2.0]),
            counts=np.array([10.0, 20.0, 30.0]),
            errors=np.array([3.16, 4.47, 5.48]),
        )
        assert spectrum.n_channels == 3
        np.testing.assert_array_equal(spectrum.energies, [0.0, 1.0, 2.0])

    def test_automatic_poisson_errors(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0, 2.0]),
            counts=np.array([10.0, 20.0, 30.0]),
        )
        expected = np.sqrt([10.0, 20.0, 30.0])
        np.testing.assert_allclose(spectrum.errors, expected)

    def test_zero_counts_get_error_one(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0]),
            counts=np.array([0.0, 100.0]),
        )
        assert spectrum.errors[0] == 1.0  # sqrt(1) for zero counts

    def test_shape_mismatch_raises(self):
        with pytest.raises(AssertionError):
            ExpSpectrum(
                energies=np.array([0.0, 1.0]),
                counts=np.array([10.0, 20.0, 30.0]),
            )

    def test_negative_energies_raises(self):
        with pytest.raises(AssertionError):
            ExpSpectrum(
                energies=np.array([-1.0, 0.0]),
                counts=np.array([10.0, 20.0]),
            )

    def test_negative_counts_raises(self):
        with pytest.raises(AssertionError):
            ExpSpectrum(
                energies=np.array([0.0, 1.0]),
                counts=np.array([10.0, -5.0]),
            )

    def test_total_counts(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0, 2.0]),
            counts=np.array([10.0, 20.0, 30.0]),
        )
        assert spectrum.total_counts == 60.0

    def test_energy_range(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.5, 1.5, 2.5]),
            counts=np.array([10.0, 20.0, 30.0]),
        )
        assert spectrum.energy_range == (0.5, 2.5)

    def test_normalize(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0, 2.0]),
            counts=np.array([10.0, 20.0, 30.0]),
        )
        norm = spectrum.normalize()
        integral = np.trapezoid(norm, spectrum.energies)
        assert np.isclose(integral, 1.0, rtol=1e-10)


class TestExpSpectrumDeadTime:
    """Dead-time correction on ExpSpectrum."""

    def test_no_correction_when_dead_time_zero(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0]),
            counts=np.array([10.0, 20.0]),
        )
        corrected = spectrum.apply_dead_time_correction()
        np.testing.assert_array_equal(corrected.counts, spectrum.counts)

    def test_correction_applied(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0]),
            counts=np.array([10.0, 20.0]),
            dead_time=0.1,
        )
        corrected = spectrum.apply_dead_time_correction()
        expected = spectrum.counts / 0.9
        np.testing.assert_allclose(corrected.counts, expected, rtol=1e-10)
        assert corrected.dead_time == 0.0

    def test_live_time_corrected_rate(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0]),
            counts=np.array([10.0, 20.0]),
            dead_time=0.1,
            live_time=10.0,
        )
        rate = spectrum.live_time_corrected_rate
        assert rate is not None
        # Total counts / (live_time * (1 - dead_time))
        expected = 30.0 / (10.0 * 0.9)
        assert np.isclose(rate, expected)


class TestExpSpectrumSerialization:
    """Dict round-trip for ExpSpectrum."""

    def test_to_from_dict(self):
        spectrum = ExpSpectrum(
            energies=np.array([0.0, 1.0, 2.0]),
            counts=np.array([10.0, 20.0, 30.0]),
            dead_time=0.05,
            live_time=100.0,
            source="Tc99",
            run_id="run_042",
            date="2025-01-15",
            metadata={"detector": "MMC", "temperature": "25 mK"},
        )

        recovered = ExpSpectrum.from_dict(spectrum.to_dict())

        np.testing.assert_allclose(recovered.energies, spectrum.energies)
        np.testing.assert_allclose(recovered.counts, spectrum.counts)
        assert recovered.dead_time == spectrum.dead_time
        assert recovered.live_time == spectrum.live_time
        assert recovered.source == spectrum.source
        assert recovered.run_id == spectrum.run_id
        assert recovered.date == spectrum.date
        assert recovered.metadata == spectrum.metadata
