"""Tests for fitter.fit_engine."""

import numpy as np

from fitter.fit_engine import SpectrumFitter, FitConfig
from fitter.model import SpectrumModel


class TestFitConfig:
    """Test FitConfig defaults."""

    def test_defaults(self):
        config = FitConfig()
        assert config.method == "trf"
        assert config.max_nfev == 1000
        assert config.loss == "linear"


class TestSpectrumFitter:
    """Test SpectrumFitter."""

    def _make_fitter(self):
        """Create a simple fitter with synthetic data."""

        # Create a simple linear model
        def linear_model(energies_keV, norm, bg):
            return norm * (1.0 - energies_keV / 100.0) + bg

        # Synthetic data
        energies = np.linspace(0.1, 99.9, 200)
        counts = 10.0 * (1.0 - energies / 100.0) + 2.0
        counts += np.random.RandomState(42).normal(0, 0.5, len(energies))

        # Mock model
        class MockModel:
            def evaluate(self, energies_keV, norm, bg):
                return linear_model(energies_keV, norm, bg)

        model = MockModel()
        fitter = SpectrumFitter(
            model, energies, counts, exp_errors=np.ones_like(counts) * 0.5
        )
        return fitter

    def test_fit_basic(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0], param_names=["norm", "bg"])

        assert result.success
        assert result.n_points == 200
        assert result.n_free == 2
        assert result.chi2 > 0
        assert result.chi2_per_dof > 0
        assert len(result.parameters) == 2

    def test_fit_with_bounds(self):
        fitter = self._make_fitter()
        result = fitter.fit(
            x0=[10.0, 2.0],
            bounds=([1.0, 0.0], [20.0, 10.0]),
        )
        assert result.success

    def test_fit_returns_residuals(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])

        assert len(result.residuals) == 200
        assert len(result.model_values) == 200

    def test_parameters_property(self):
        fitter = self._make_fitter()
        fitter.fit(x0=[10.0, 2.0])
        params = fitter.parameters
        assert len(params) == 2

    def test_covariance_property(self):
        fitter = self._make_fitter()
        fitter.fit(x0=[10.0, 2.0])
        cov = fitter.covariance
        assert cov.shape == (2, 2)

    def test_chi2_per_dof(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        assert result.chi2_per_dof == result.chi2 / result.n_free

    def test_p_value(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        assert 0.0 <= result.p_value <= 1.0

    def test_parameters_with_errors(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        params = result.parameters_with_errors
        assert len(params) == 2
        for name, (val, err) in params.items():
            assert isinstance(val, float)
            assert isinstance(err, float)
            assert err >= 0

    def test_correlation_matrix(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        corr = result.correlation_matrix
        assert corr.shape == (2, 2)
        # Diagonal should be 1
        np.testing.assert_allclose(np.diag(corr), 1.0)

    def test_summary(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0], param_names=["norm", "bg"])
        summary = result.summary(param_names=["norm", "bg"])
        assert "norm" in summary
        assert "bg" in summary
        assert "Fit Results" in summary

    def test_profile_likelihood(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        param_vals, chi2_vals, norm_vals = result.profile_likelihood(0, n_points=10)
        assert len(param_vals) == 10
        assert len(chi2_vals) == 10
        assert len(norm_vals) == 10
        assert np.min(norm_vals) >= 0  # Normalized χ² should be >= 0

    def test_confidence_interval(self):
        fitter = self._make_fitter()
        result = fitter.fit(x0=[10.0, 2.0])
        ci = result.confidence_interval(0, confidence=0.95)
        assert len(ci) == 2
        assert ci[0] <= ci[1]


class TestSpectrumFitterWithRealModel:
    """Test SpectrumFitter with a real SpectrumModel."""

    def test_fit_with_real_model(self):
        from beta_spectrum import BetaSpectrum, SpectrumConfig

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            detector_model="gaussian",
            detector_sigma_a_keV=0.5,
            detector_n_channels=200,
        )
        spectrum = BetaSpectrum.from_config(config)
        detector = BetaSpectrum.create_detector_from_config(config)
        model = SpectrumModel(spectrum, detector)

        # Generate synthetic data from model
        energies_keV = detector.channel_energies
        true_norm = 5.0
        true_bg = 0.1
        counts = model.evaluate(energies_keV, true_norm, true_bg)
        counts = np.maximum(counts, 0)  # Clip negatives
        counts += np.random.RandomState(42).normal(
            0, np.sqrt(np.maximum(counts, 1)), len(energies_keV)
        )
        counts = np.maximum(counts, 0)

        fitter = SpectrumFitter(
            model,
            energies_keV,
            counts,
            exp_errors=np.sqrt(np.maximum(counts, 1)),
        )
        result = fitter.fit(x0=[1.0, 0.0], param_names=["norm", "bg"])

        assert result.success
        # Fit should complete without error
        assert len(result.parameters) == 2
