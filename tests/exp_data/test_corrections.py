"""Tests for exp_data.corrections."""

import numpy as np

from exp_data.corrections import (
    DeadTimeCorrection,
    PileUpCorrection,
    BackgroundSubtractor,
)


class TestDeadTimeCorrection:
    """Test dead-time correction."""

    def test_non_paralyzable_live_time(self):
        dt = DeadTimeCorrection(dead_time=0.001, model="non_paralyzable")
        # At 100 cps: live = 1 - 100*0.001 = 0.9
        assert np.isclose(dt.live_time_fraction(100.0), 0.9)

    def test_non_paralyzable_correct_rate(self):
        dt = DeadTimeCorrection(dead_time=0.01, model="non_paralyzable")
        # True rate = observed / (1 - obs*dt) = 90 / (1 - 90*0.01) = 90 / 0.1 = 900
        assert np.isclose(dt.correct_rate(90.0), 900.0, rtol=0.01)

    def test_paralyzable_live_time(self):
        dt = DeadTimeCorrection(dead_time=0.001, model="paralyzable")
        # At low rate: exp(-0.1) ≈ 0.905
        frac = dt.live_time_fraction(100.0)
        assert np.isclose(frac, np.exp(-0.1), rtol=1e-5)

    def test_correct_counts(self):
        dt = DeadTimeCorrection(dead_time=0.01, model="non_paralyzable")
        counts = np.array([100.0, 200.0, 300.0])
        live_time = 10.0
        corrected = dt.correct_counts(counts, live_time)
        # Total observed rate = 600/10 = 60 cps
        # True rate = 60 / (1 - 60*0.01) = 60 / 0.94 = 63.83 cps
        # Scale factor = 63.83/60 = 1.0638
        scale = 1.0 / (1 - 60.0 * 0.01)
        expected = counts * scale
        np.testing.assert_allclose(corrected, expected)


class TestPileUpCorrection:
    """Test pile-up correction."""

    def test_pile_up_fraction(self):
        pu = PileUpCorrection(dead_time=0.001)
        # At 100 cps: (0.1)² = 0.01
        frac = pu.estimate_pile_up_fraction(100.0)
        assert np.isclose(frac, 0.01)

    def test_correct_spectrum(self):
        pu = PileUpCorrection(dead_time=0.001)
        counts = np.array([100.0, 200.0, 300.0])
        live_time = 10.0
        corrected = pu.correct_spectrum(counts, live_time)
        # Pile-up correction adds a small fraction back
        assert np.all(corrected >= counts)


class TestBackgroundSubtractor:
    """Test background subtraction."""

    def test_constant_background(self):
        energies = np.linspace(0, 100, 200)
        counts = 50 * np.exp(-0.5 * ((energies - 50) / 5) ** 2) + 10.0  # peak + bg

        corrected = BackgroundSubtractor.constant_background(
            counts,
            energies,
            low_roi=(0.0, 10.0),
            high_roi=(90.0, 100.0),
        )

        # Peak should still be visible, background removed
        assert np.max(corrected) > 0
        # Background regions should be near zero
        assert np.mean(corrected[energies <= 10]) < 5
        assert np.mean(corrected[energies >= 90]) < 5

    def test_polynomial_background(self):
        energies = np.linspace(0, 100, 200)
        counts = 50 * np.exp(-0.5 * ((energies - 50) / 5) ** 2) + 10.0 + 0.05 * energies

        corrected = BackgroundSubtractor.polynomial_background(
            counts,
            energies,
            order=1,
            low_roi=(0.0, 10.0),
            high_roi=(90.0, 100.0),
        )

        assert np.max(corrected) > 0

    def test_roi_average(self):
        counts = np.array([10.0, 12.0, 50.0, 10.0, 11.0])
        mask = np.array([True, True, False, True, True])
        avg = BackgroundSubtractor.roi_average(counts, mask)
        np.testing.assert_allclose(avg, 10.75)

    def test_roi_average_empty(self):
        counts = np.array([10.0, 20.0])
        mask = np.array([False, False])
        avg = BackgroundSubtractor.roi_average(counts, mask)
        assert avg == 0.0
