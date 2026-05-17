"""Tests for exp_data.fitters."""

import numpy as np

from exp_data.fitters import GaussianFitter, PeakFitter


class TestGaussianFitter:
    """Test Gaussian peak fitting."""

    def _make_synthetic_peak(
        self,
        center=50.0,
        sigma=2.0,
        amplitude=100.0,
        background=10.0,
        n_points=500,
        energy_range=(0, 100),
    ):
        """Create synthetic Gaussian peak with noise."""
        energies = np.linspace(energy_range[0], energy_range[1], n_points)
        counts = (
            amplitude * np.exp(-0.5 * ((energies - center) / sigma) ** 2) + background
        )
        counts += np.random.RandomState(42).normal(
            0, np.sqrt(np.maximum(counts, 1)), n_points
        )
        return energies, counts

    def test_fit_returns_valid_result(self):
        energies, counts = self._make_synthetic_peak()
        fitter = GaussianFitter(energies, counts)
        result = fitter.fit()

        assert result.success
        assert result.chi2 > 0
        assert result.n_dof > 0
        assert result.covariance.shape == (4, 4)

    def test_fit_recovers_parameters(self):
        center, sigma, amplitude, background = 50.0, 2.0, 100.0, 10.0
        energies, counts = self._make_synthetic_peak(
            center=center, sigma=sigma, amplitude=amplitude, background=background
        )
        fitter = GaussianFitter(energies, counts)
        result = fitter.fit()

        np.testing.assert_allclose(result.center_keV, center, rtol=0.05)
        np.testing.assert_allclose(result.sigma_keV, sigma, rtol=0.1)
        np.testing.assert_allclose(result.amplitude, amplitude, rtol=0.15)

    def test_resolution_property(self):
        energies, counts = self._make_synthetic_peak(sigma=2.0)
        fitter = GaussianFitter(energies, counts)
        result = fitter.fit()

        expected_fwhm = 2.35482 * 2.0
        np.testing.assert_allclose(result.resolution(), expected_fwhm, rtol=0.05)

    def test_resolution_fraction(self):
        energies, counts = self._make_synthetic_peak(center=50.0, sigma=2.0)
        fitter = GaussianFitter(energies, counts)
        result = fitter.fit()

        expected_frac = (2.35482 * 2.0) / 50.0
        np.testing.assert_allclose(
            result.resolution_fraction(), expected_frac, rtol=0.05
        )

    def test_auto_guess(self):
        energies, counts = self._make_synthetic_peak()
        fitter = GaussianFitter(energies, counts)
        result = fitter.fit()
        assert result.success


class TestPeakFitter:
    """Test multi-peak detection and fitting."""

    def _make_multi_peak(self):
        """Create spectrum with two peaks."""
        energies = np.linspace(0, 100, 1000)
        counts = np.zeros_like(energies)
        counts += 50 * np.exp(-0.5 * ((energies - 20) / 2) ** 2)
        counts += 80 * np.exp(-0.5 * ((energies - 70) / 3) ** 2)
        counts += 5  # constant background
        counts += np.random.RandomState(42).normal(
            0, np.sqrt(np.maximum(counts, 1)), len(energies)
        )
        return energies, counts

    def test_find_peaks(self):
        energies, counts = self._make_multi_peak()
        fitter = PeakFitter(energies, counts)
        peaks = fitter.find_peaks(min_height=10.0, min_distance_keV=10.0)

        assert len(peaks) >= 1  # At least one peak should be found

    def test_fit_peaks_with_known_energies(self):
        energies, counts = self._make_multi_peak()
        fitter = PeakFitter(energies, counts)
        results = fitter.fit_peaks(known_energies=np.array([20.0, 70.0]))

        assert len(results) >= 1
        for r in results:
            assert r.success
            assert r.n_dof > 0

    def test_fit_peaks_auto_detect(self):
        energies, counts = self._make_multi_peak()
        fitter = PeakFitter(energies, counts)
        results = fitter.fit_peaks()

        assert len(results) >= 1
