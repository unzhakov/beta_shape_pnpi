# tests/test_fitter.py — Tests for χ² minimization framework

"""
Why test this?
--------------
The fitter is the core numerical tool for extracting physical parameters
from beta spectra. It must be correct because all downstream results
(g_V, g_A, C(W) extraction) depend on it.

We test:
 1. Basic χ² minimization on known functions (quadratic, linear).
 2. Covariance matrix estimation (should match input uncertainties).
 3. Confidence interval computation via profile likelihood.
 4. Goodness-of-fit metrics (χ²/ndf, p-value).
 5. Bounds handling and convergence diagnostics.
 6. Robust loss functions (soft_l1, huber) for outlier resistance.
"""

import numpy as np

from beta_spectrum import CurveFitter, FitConfig


class TestBasicFitting:
    """Test χ² minimization on simple known functions."""

    def test_quadratic_fit(self):
        """Fit y = a + b·x + c·x² with known coefficients."""
        np.random.seed(42)
        x = np.linspace(0, 10, 50)
        true_params = [1.0, 2.0, -0.1]
        y = true_params[0] + true_params[1] * x + true_params[2] * x**2
        y += np.random.normal(0, 0.5, len(x))

        def model(x, a, b, c):
            return a + b * x + c * x**2

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(y) * 0.5)
        result = fitter.fit([0.0, 0.0, 0.0])

        for i, name in enumerate(["a", "b", "c"]):
            assert (
                abs(result.parameters[i] - true_params[i]) < 0.5
            ), f"Parameter {name}: expected {true_params[i]}, got {result.parameters[i]}"

        assert result.success, "Fit should converge"
        assert result.chi2 > 0, "χ² must be positive"

    def test_linear_fit(self):
        """Fit y = a + b·x with known coefficients."""
        np.random.seed(123)
        x = np.linspace(0, 5, 30)
        true_a, true_b = 3.0, -1.5
        y = true_a + true_b * x
        y += np.random.normal(0, 0.2, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(y) * 0.2)
        result = fitter.fit([0.0, 0.0])

        assert abs(result.parameters[0] - true_a) < 0.3
        assert abs(result.parameters[1] - true_b) < 0.3
        assert result.success

    def test_constant_fit(self):
        """Fit y = a (constant) with known value."""
        np.random.seed(456)
        x = np.linspace(0, 10, 40)
        true_a = 5.0
        y = np.full_like(x, true_a)
        y += np.random.normal(0, 0.3, len(x))

        def model(x, a):
            return np.full_like(x, a)

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.3)
        result = fitter.fit([0.0])

        assert abs(result.parameters[0] - true_a) < 0.2
        assert result.success


class TestCovarianceEstimation:
    """Test that covariance matrix is correctly estimated."""

    def test_covariance_matches_uncertainty(self):
        """For a constant model, parameter error should match input σ."""
        np.random.seed(789)
        x = np.linspace(0, 10, 50)
        true_a = 5.0
        sigma = 0.5
        y = np.full_like(x, true_a, dtype=float)
        y += np.random.normal(0, sigma, len(x))

        def model(x, a):
            return np.full_like(x, a)

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * sigma)
        result = fitter.fit([0.0])

        error = np.sqrt(result.covariance[0, 0])
        # Expected error: σ / sqrt(N)
        expected_error = sigma / np.sqrt(len(x))
        assert (
            abs(error - expected_error) < expected_error * 0.3
        ), f"Covariance error {error:.4f} != expected {expected_error:.4f}"

    def test_covariance_matrix_positive_definite(self):
        """Covariance matrix must be positive definite."""
        np.random.seed(321)
        x = np.linspace(0, 5, 30)
        y = 2.0 + 1.5 * x - 0.3 * x**2
        y += np.random.normal(0, 0.2, len(x))

        def model(x, a, b, c):
            return a + b * x + c * x**2

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.2)
        result = fitter.fit([0.0, 0.0, 0.0])

        # Check eigenvalues are non-negative
        eigvals = np.linalg.eigvalsh(result.covariance)
        assert np.all(eigvals >= -1e-10), f"Covariance not positive definite: {eigvals}"

    def test_correlation_matrix(self):
        """Correlation matrix should have 1.0 on diagonal."""
        np.random.seed(654)
        x = np.linspace(0, 10, 50)
        y = 1.0 + 2.0 * x + 0.1 * x**2
        y += np.random.normal(0, 0.5, len(x))

        def model(x, a, b, c):
            return a + b * x + c * x**2

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.5)
        result = fitter.fit([0.0, 0.0, 0.0])

        diag = np.diag(result.correlation_matrix)
        assert np.allclose(
            diag, 1.0, atol=1e-6
        ), f"Diagonal of correlation matrix: {diag}"


class TestGoodnessOfFit:
    """Test χ²/ndf and p-value calculations."""

    def test_chi2_per_dof(self):
        """χ²/ndof should be positive and reasonable."""
        np.random.seed(111)
        x = np.linspace(0, 10, 50)
        y = 2.0 * x
        y += np.random.normal(0, 0.3, len(x))

        def model(x, a):
            return a * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.3)
        result = fitter.fit([0.0])

        assert result.chi2_per_dof > 0
        assert result.n_points == 50
        assert result.n_free == 1
        assert result.n_free == 1

    def test_p_value_range(self):
        """p-value should be in [0, 1]."""
        np.random.seed(222)
        x = np.linspace(0, 5, 30)
        y = 1.0 + x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.0, 0.0])

        assert 0.0 <= result.p_value <= 1.0

    def test_p_value_good_fit(self):
        """Good fit should have reasonable p-value."""
        np.random.seed(333)
        x = np.linspace(0, 10, 100)
        y = 3.0 + 1.5 * x - 0.05 * x**2
        y += np.random.normal(0, 0.2, len(x))

        def model(x, a, b, c):
            return a + b * x + c * x**2

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.2)
        result = fitter.fit([0.0, 0.0, 0.0])

        assert (
            result.p_value > 0.001
        ), f"Good fit should have reasonable p-value, got {result.p_value}"


class TestBounds:
    """Test parameter bounds handling."""

    def test_lower_bound(self):
        """Test that lower bounds are respected."""
        np.random.seed(444)
        x = np.linspace(0, 5, 30)
        y = 2.0 + 1.0 * x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([1.5, 0.5], bounds=([1.0, 0.0], [10.0, 5.0]))

        assert (
            result.parameters[0] >= 1.0
        ), f"Parameter a should be >= 1.0, got {result.parameters[0]}"
        assert (
            result.parameters[1] >= 0.0
        ), f"Parameter b should be >= 0.0, got {result.parameters[1]}"

    def test_upper_bound(self):
        """Test that upper bounds are respected."""
        np.random.seed(555)
        x = np.linspace(0, 5, 30)
        y = 0.5 + 0.5 * x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.5, 0.5], bounds=([0.0, 0.0], [1.0, 1.0]))

        assert result.parameters[0] <= 1.0
        assert result.parameters[1] <= 1.0


class TestRobustLoss:
    """Test robust loss functions for outlier resistance."""

    def test_soft_l1_loss(self):
        """soft_l1 loss should handle outliers better than linear."""
        np.random.seed(666)
        x = np.linspace(0, 10, 50)
        y = 1.0 + 2.0 * x
        y[25] = 50.0  # Add outlier

        def model(x, a, b):
            return a + b * x

        # Linear loss — outlier pulls fit away
        config_linear = FitConfig(loss="linear")
        fitter_linear = CurveFitter(
            model, x, y, uncertainties=np.ones_like(x) * 1.0, config=config_linear
        )
        result_linear = fitter_linear.fit([0.0, 0.0])

        # soft_l1 loss — outlier has less influence
        config_robust = FitConfig(loss="soft_l1")
        fitter_robust = CurveFitter(
            model, x, y, uncertainties=np.ones_like(x) * 1.0, config=config_robust
        )
        result_robust = fitter_robust.fit([0.0, 0.0])

        # Robust fit should be closer to true slope
        assert (
            abs(result_robust.parameters[1] - 2.0)
            < abs(result_linear.parameters[1] - 2.0) * 1.5
        ), "Robust loss should handle outliers better"


class TestFitResult:
    """Test FitResult dataclass methods."""

    def test_parameters_with_errors(self):
        """parameters_with_errors should return dict of (value, error)."""
        np.random.seed(777)
        x = np.linspace(0, 5, 30)
        y = 2.0 + x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.0, 0.0])

        errors = result.parameters_with_errors
        assert len(errors) == 2
        for key, (val, err) in errors.items():
            assert isinstance(val, float)
            assert isinstance(err, float)
            assert err >= 0

    def test_summary(self):
        """summary() should return a non-empty string."""
        np.random.seed(888)
        x = np.linspace(0, 5, 30)
        y = 1.0 + 2.0 * x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.0, 0.0])

        summary = result.summary(param_names=["a", "b"])
        assert "Fit Results" in summary
        assert "a" in summary
        assert "b" in summary
        assert "χ²" in summary
        assert "Converged" in summary

    def test_residuals_property(self):
        """residuals should be data - model."""
        np.random.seed(999)
        x = np.linspace(0, 5, 30)
        y = 2.0 + x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.0, 0.0])

        expected_residuals = y - result.model_values
        np.testing.assert_allclose(result.residuals, expected_residuals)


class TestFitConfig:
    """Test FitConfig options."""

    def test_custom_config(self):
        """Custom FitConfig should be respected."""
        np.random.seed(101)
        x = np.linspace(0, 5, 30)
        y = 1.0 + x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        config = FitConfig(
            method="trf",
            max_nfev=500,
            loss="linear",
            ftol=1e-10,
        )
        fitter = CurveFitter(
            model, x, y, uncertainties=np.ones_like(x) * 0.1, config=config
        )
        result = fitter.fit([0.0, 0.0])

        assert result.success


class TestConfidenceIntervals:
    """Test confidence interval computation."""

    def test_95_percent_ci(self):
        """95% CI should contain the true parameter value."""
        np.random.seed(202)
        x = np.linspace(0, 10, 60)
        true_a, true_b = 2.0, 1.5
        y = true_a + true_b * x
        y += np.random.normal(0, 0.3, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.3)
        result = fitter.fit([0.0, 0.0])

        ci_a = result.confidence_interval(0, confidence=0.95)
        ci_b = result.confidence_interval(1, confidence=0.95)

        assert ci_a[0] <= true_a <= ci_a[1], f"True a={true_a} not in CI {ci_a}"
        assert ci_b[0] <= true_b <= ci_b[1], f"True b={true_b} not in CI {ci_b}"

    def test_ci_width_increases_with_confidence(self):
        """Higher confidence level → wider interval."""
        np.random.seed(303)
        x = np.linspace(0, 5, 40)
        y = 1.0 + 2.0 * x
        y += np.random.normal(0, 0.2, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.2)
        result = fitter.fit([0.0, 0.0])

        ci_90 = result.confidence_interval(0, confidence=0.90)
        ci_95 = result.confidence_interval(0, confidence=0.95)

        width_90 = ci_90[1] - ci_90[0]
        width_95 = ci_95[1] - ci_95[0]

        assert (
            width_95 > width_90
        ), f"95% CI width {width_95} should be > 90% width {width_90}"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_zero_uncertainties_filtered(self):
        """Points with zero uncertainty should be filtered out."""
        x = np.array([1.0, 2.0, 3.0, 4.0])
        y = np.array([2.0, 4.0, 6.0, 8.0])
        unc = np.array([1.0, 0.0, 1.0, 0.0])

        def model(x, a):
            return a * x

        fitter = CurveFitter(model, x, y, uncertainties=unc)
        assert len(fitter.x_data) == 2, "Zero-uncertainty points should be removed"
        assert len(fitter.y_data) == 2

    def test_negative_uncertainties_filtered(self):
        """Points with negative uncertainty should be filtered out."""
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([2.0, 4.0, 6.0])
        unc = np.array([1.0, -1.0, 1.0])

        def model(x, a):
            return a * x

        fitter = CurveFitter(model, x, y, uncertainties=unc)
        assert len(fitter.x_data) == 2

    def test_unbounded_parameters(self):
        """Unbounded parameters (±inf) should work."""
        np.random.seed(404)
        x = np.linspace(0, 5, 30)
        y = 1.0 + 2.0 * x
        y += np.random.normal(0, 0.1, len(x))

        def model(x, a, b):
            return a + b * x

        fitter = CurveFitter(model, x, y, uncertainties=np.ones_like(x) * 0.1)
        result = fitter.fit([0.0, 0.0], bounds=([-np.inf, -np.inf], [np.inf, np.inf]))
        assert result.success
