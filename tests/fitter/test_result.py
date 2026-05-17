"""Tests for fitter.result."""

import numpy as np

from fitter.result import AnalysisResult, FitSummary


class MockFitResult:
    """Mock FitResult for testing AnalysisResult."""

    def __init__(self):
        self.success = True
        self.chi2 = 12.5
        self.n_points = 100
        self.n_free = 3
        self.parameters = np.array([5.0, 0.1, 2.0])
        self.covariance = np.array(
            [
                [0.01, 0.0, 0.0],
                [0.0, 0.001, 0.0],
                [0.0, 0.0, 0.04],
            ]
        )
        self.model_values = np.ones(100)
        self.residuals = np.zeros(100)
        self.message = "Optimization terminated successfully."

    @property
    def chi2_per_dof(self):
        return self.chi2 / (self.n_points - self.n_free)

    @property
    def p_value(self):
        from scipy.stats import chi2

        return 1.0 - chi2.cdf(self.chi2, self.n_points - self.n_free)

    @property
    def parameters_with_errors(self):
        errors = np.sqrt(np.diag(self.covariance))
        return {
            f"p{i}": (float(self.parameters[i]), float(errors[i]))
            for i in range(len(self.parameters))
        }


class TestFitSummary:
    """Test FitSummary."""

    def test_str(self):
        summary = FitSummary(
            nuclide="Tc99",
            converged=True,
            chi2_per_dof=1.234,
            p_value=0.0567,
            parameters={"norm": (5.0, 0.1), "bg": (0.1, 0.01)},
        )
        text = str(summary)
        assert "Tc99" in text
        assert "1.234" in text
        assert "5.0" in text
        assert "Converged: True" in text


class TestAnalysisResult:
    """Test AnalysisResult."""

    def test_summary(self):
        mock = MockFitResult()
        result = AnalysisResult(
            fit_result=mock,
            nuclide="Tc99",
            metadata={"param_names": ["norm", "bg", "endpoint"]},
        )
        summary = result.summary
        assert summary.nuclide == "Tc99"
        assert summary.converged is True
        assert summary.chi2_per_dof > 0

    def test_to_dict(self):
        mock = MockFitResult()
        result = AnalysisResult(
            fit_result=mock,
            nuclide="Tc99",
            gv_eff=0.97,
            ga_eff=0.57,
            endpoint_keV=294.0,
        )
        d = result.to_dict()
        assert d["nuclide"] == "Tc99"
        assert d["gv_eff"] == 0.97
        assert d["ga_eff"] == 0.57
        assert d["endpoint_keV"] == 294.0
        assert d["converged"] is True

    def test_empty_result(self):
        mock = MockFitResult()
        result = AnalysisResult(fit_result=mock)
        summary = result.summary
        assert summary.nuclide == ""
