"""
Multi-parameter fitting engine for beta spectrum analysis.

Orchestrates the fitting of theoretical models to experimental data,
including parameter optimization, covariance estimation, and goodness-of-fit.

Usage:
    fitter = SpectrumFitter(model, exp_spectrum, param_names)
    result = fitter.fit(x0=[norm, bg], bounds=([0, 0], [100, 1e6]))
    print(result.summary())
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares


@dataclass
class FitConfig:
    """Configuration for spectrum fitting."""

    method: str = "trf"
    """Optimization method: 'trf', 'dogbox', 'lm'."""

    max_nfev: int = 1000
    """Maximum function evaluations."""

    x_scale: str | float = "jac"
    """Parameter scaling."""

    loss: str = "linear"
    """Loss function: 'linear', 'soft_l1', 'huber', 'cauchy'."""

    f_scale: float = 1.0
    """Robust loss scale parameter."""

    ftol: float = 1e-10
    """Function value tolerance."""

    xtol: float = 1e-10
    """Parameter change tolerance."""

    gtol: float = 1e-10
    """Gradient tolerance."""


@dataclass
class FitResult:
    """Container for fit results."""

    parameters: np.ndarray
    """Best-fit parameter values."""

    covariance: np.ndarray
    """Parameter covariance matrix."""

    chi2: float
    """Minimum χ² value."""

    n_points: int
    """Number of data points."""

    n_free: int
    """Number of free parameters."""

    model_values: np.ndarray
    """Model values at best-fit parameters."""

    residuals: np.ndarray
    """Data minus model."""

    success: bool
    """Whether optimization converged."""

    message: str
    """Optimization message."""

    _residual_fn: Optional[Callable] = field(default=None, repr=False)
    """Residual function for profile likelihood."""

    _x_data: Optional[np.ndarray] = field(default=None, repr=False)
    """Independent variable data."""

    _y_data: Optional[np.ndarray] = field(default=None, repr=False)
    """Dependent variable data."""

    _uncertainties: Optional[np.ndarray] = field(default=None, repr=False)
    """Measurement uncertainties."""

    @property
    def chi2_per_dof(self) -> float:
        """χ² per degree of freedom."""
        return self.chi2 / max(self.n_free, 1)

    @property
    def p_value(self) -> float:
        """p-value for the fit (survival function of χ²)."""
        from scipy.stats import chi2

        return float(1.0 - chi2.cdf(self.chi2, self.n_points - self.n_free))

    @property
    def parameters_with_errors(self) -> dict[str, tuple[float, float]]:
        """Parameter values with 1σ uncertainties."""
        errors = np.sqrt(np.maximum(np.diag(self.covariance), 0.0))
        return {
            f"p{i}": (float(self.parameters[i]), float(errors[i]))
            for i in range(len(self.parameters))
        }

    @property
    def correlation_matrix(self) -> np.ndarray:
        """Parameter correlation matrix."""
        diag = np.sqrt(np.maximum(np.diag(self.covariance), 0.0))
        if np.any(diag == 0):
            return np.zeros_like(self.covariance)
        return self.covariance / np.outer(diag, diag)

    def summary(self, param_names: Optional[List[str]] = None) -> str:
        """Human-readable fit summary."""
        lines = [
            "=" * 60,
            "Fit Results",
            "=" * 60,
            f"Converged: {self.success}",
            f"Message: {self.message}",
            f"χ² = {self.chi2:.4f}",
            f"χ²/ndof = {self.chi2_per_dof:.3f} (ndof = {self.n_free})",
            f"p-value = {self.p_value:.4f}",
            "-" * 60,
            "Parameters:",
        ]

        errors = np.sqrt(np.maximum(np.diag(self.covariance), 0.0))
        for i in range(len(self.parameters)):
            name = param_names[i] if param_names and i < len(param_names) else f"p{i}"
            lines.append(
                f"  {name:20s} = {self.parameters[i]:14.6f}  ±  {errors[i]:10.6f}"
            )

        lines.append("-" * 60)
        lines.append("Correlation matrix (off-diagonal):")
        corr = self.correlation_matrix
        n = len(self.parameters)
        header = "          " + "".join(f"{f'p{i}':>12s}" for i in range(n))
        lines.append(header)
        for i in range(n):
            name = param_names[i] if param_names and i < len(param_names) else f"p{i}"
            row = f"  {name:8s}" + "".join(
                f"{corr[i, j]:12.3f}" for j in range(n) if j != i
            )
            lines.append(row)

        lines.append("=" * 60)
        return "\n".join(lines)

    def profile_likelihood(
        self,
        param_index: int,
        n_points: int = 50,
        param_range: Optional[Tuple[float, float]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute profile likelihood for a single parameter."""
        if self._residual_fn is None:
            raise RuntimeError("Profile likelihood not available.")

        if param_range is None:
            sigma = np.sqrt(abs(self.covariance[param_index, param_index]))
            center = self.parameters[param_index]
            param_range = (center - 3 * sigma, center + 3 * sigma)

        param_values = np.linspace(param_range[0], param_range[1], n_points)
        chi2_values = np.zeros(n_points)

        for i, pval in enumerate(param_values):
            params = list(self.parameters)
            params[param_index] = pval
            res = self._residual_fn(np.array(params))
            chi2_values[i] = np.sum(res**2)

        min_chi2 = np.min(chi2_values)
        normalized = chi2_values - min_chi2

        return param_values, chi2_values, normalized

    def confidence_interval(
        self,
        param_index: int,
        confidence: float = 0.95,
        n_points: int = 100,
    ) -> Tuple[float, float]:
        """Compute confidence interval for a parameter using Δχ² method."""
        from scipy.stats import chi2

        delta_chi2 = chi2.ppf(confidence, df=1)
        param_values, chi2_values, normalized = self.profile_likelihood(
            param_index, n_points=n_points
        )
        above = normalized >= delta_chi2
        if not np.any(above):
            return (
                float(self.parameters[param_index]),
                float(self.parameters[param_index]),
            )
        crossings = np.where(np.diff(above.astype(int)))[0]
        if len(crossings) >= 2:
            lower = param_values[crossings[0]]
            upper = param_values[crossings[-1] + 1]
        elif len(crossings) == 1:
            lower = param_values[crossings[0]]
            upper = param_values[-1]
        else:
            lower = param_values[0]
            upper = param_values[-1]
        return (float(lower), float(upper))


class SpectrumFitter:
    """Multi-parameter fitting engine for beta spectrum analysis.

    Wraps a :class:`SpectrumModel` and experimental data to perform
    χ² minimization with configurable loss functions and bounds.

    Parameters
    ----------
    model : SpectrumModel
        Theoretical model (with optional detector convolution).
    exp_energies : np.ndarray
        Experimental energy grid (keV).
    exp_counts : np.ndarray
        Experimental counts per channel.
    exp_errors : np.ndarray, optional
        Experimental uncertainties. If None, assumed to be √counts.
    config : FitConfig, optional
        Fitting configuration.
    """

    def __init__(
        self,
        model: Any,
        exp_energies: np.ndarray,
        exp_counts: np.ndarray,
        exp_errors: Optional[np.ndarray] = None,
        config: Optional[FitConfig] = None,
    ):
        self.model = model
        self.exp_energies = np.asarray(exp_energies, dtype=np.float64)
        self.exp_counts = np.asarray(exp_counts, dtype=np.float64)

        if exp_errors is not None:
            self.exp_errors = np.asarray(exp_errors, dtype=np.float64)
        else:
            self.exp_errors = np.sqrt(np.maximum(self.exp_counts, 1.0))

        # Filter out zero-uncertainty points
        valid = self.exp_errors > 0
        self.exp_energies = self.exp_energies[valid]
        self.exp_counts = self.exp_counts[valid]
        self.exp_errors = self.exp_errors[valid]

        self.config = config or FitConfig()
        self._n_params: Optional[int] = None
        self._last_result: Optional[FitResult] = None

    def _residuals(self, params: np.ndarray) -> np.ndarray:
        """Compute weighted residuals."""
        model_values = self.model.evaluate(self.exp_energies, params[0], params[1])
        return (self.exp_counts - model_values) / self.exp_errors

    def _chi2(self, params: np.ndarray) -> float:
        """Compute χ²."""
        model_values = self.model.evaluate(self.exp_energies, params[0], params[1])
        return float(np.sum(((self.exp_counts - model_values) / self.exp_errors) ** 2))

    def fit(
        self,
        x0: Optional[Sequence[float]] = None,
        bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
        param_names: Optional[List[str]] = None,
    ) -> FitResult:
        """Perform χ² minimization.

        Parameters
        ----------
        x0 : sequence of float, optional
            Initial parameter values. Auto-estimated if None.
        bounds : (lower, upper), optional
            Parameter bounds.
        param_names : list of str, optional
            Parameter names for summary output.

        Returns
        -------
        FitResult
        """
        if x0 is None:
            x0 = [1.0, 0.0]

        x0_array = np.asarray(x0, dtype=np.float64)
        self._n_params = len(x0_array)

        bounds_array: Optional[Tuple[np.ndarray, np.ndarray]] = None
        if bounds is not None:
            bounds_array = (
                np.asarray(bounds[0], dtype=np.float64),
                np.asarray(bounds[1], dtype=np.float64),
            )

        ls_kwargs: dict[str, Any] = {
            "fun": self._residuals,
            "x0": x0_array,
            "method": self.config.method,
            "max_nfev": self.config.max_nfev,
            "x_scale": self.config.x_scale,
            "loss": self.config.loss,
            "f_scale": self.config.f_scale,
            "ftol": self.config.ftol,
            "xtol": self.config.xtol,
            "gtol": self.config.gtol,
        }
        if bounds is not None:
            ls_kwargs["bounds"] = bounds_array

        result = least_squares(**ls_kwargs)

        # Compute covariance
        cov = self._compute_covariance(result)

        model_values = self.model.evaluate(self.exp_energies, *result.x)
        residuals = self.exp_counts - model_values

        fit_result = FitResult(
            parameters=result.x,
            covariance=cov,
            chi2=self._chi2(result.x),
            n_points=len(self.exp_energies),
            n_free=self._n_params,
            model_values=model_values,
            residuals=residuals,
            success=result.success,
            message=result.message,
            _residual_fn=self._residuals,
            _x_data=self.exp_energies,
            _y_data=self.exp_counts,
            _uncertainties=self.exp_errors,
        )

        self._last_result = fit_result
        return fit_result

    def _compute_covariance(self, result: Any) -> np.ndarray:
        """Estimate parameter covariance from Jacobian."""
        n = len(result.x)
        cov = np.zeros((n, n))

        if result.jac is not None:
            jtj = result.jac.T @ result.jac
            if np.linalg.cond(jtj) < 1e12:
                cov = np.linalg.inv(jtj)
                chi2 = self._chi2(result.x)
                ndof = max(len(self.exp_energies) - n, 1)
                cov = cov * chi2 / ndof
            else:
                cov = np.linalg.pinv(jtj)
                chi2 = self._chi2(result.x)
                ndof = max(len(self.exp_energies) - n, 1)
                cov = cov * chi2 / ndof

        return cov

    @property
    def parameters(self) -> np.ndarray:
        """Current best-fit parameters."""
        if self._last_result is None:
            raise RuntimeError("No fit result. Call fit() first.")
        return self._last_result.parameters

    @property
    def covariance(self) -> np.ndarray:
        """Current covariance matrix."""
        if self._last_result is None:
            raise RuntimeError("No fit result. Call fit() first.")
        return self._last_result.covariance


class CurveFitter:
    """General-purpose χ² curve fitter for arbitrary callable models.

    This is a thin wrapper around scipy's least_squares that provides
    the same API as the original beta_spectrum.fitter.CurveFitter.

    Parameters
    ----------
    model : callable
        Function f(x, *params) that computes model values.
    x_data : np.ndarray
        Independent variable data.
    y_data : np.ndarray
        Dependent variable data.
    uncertainties : np.ndarray, optional
        Measurement uncertainties. If None, assumed to be 1.0.
    config : FitConfig, optional
        Fitting configuration.

    Examples
    --------
    >>> def linear(x, a, b):
    ...     return a + b * x
    >>> fitter = CurveFitter(linear, x, y, uncertainties=unc)
    >>> result = fitter.fit([0.0, 0.0])
    """

    def __init__(
        self,
        model: Callable,
        x_data: np.ndarray,
        y_data: np.ndarray,
        uncertainties: Optional[np.ndarray] = None,
        config: Optional[FitConfig] = None,
    ):
        self.model = model
        self.x_data = np.asarray(x_data, dtype=np.float64)
        self.y_data = np.asarray(y_data, dtype=np.float64)

        if uncertainties is not None:
            self.uncertainties = np.asarray(uncertainties, dtype=np.float64)
        else:
            self.uncertainties = np.ones_like(self.y_data)

        # Filter out zero/negative uncertainty points
        valid = self.uncertainties > 0
        self.x_data = self.x_data[valid]
        self.y_data = self.y_data[valid]
        self.uncertainties = self.uncertainties[valid]

        self.config = config or FitConfig()
        self._last_result: Optional[FitResult] = None
        self._n_params: Optional[int] = None

    def _residuals(self, params: np.ndarray) -> np.ndarray:
        """Compute weighted residuals."""
        model_values = self.model(self.x_data, *params)
        return (self.y_data - model_values) / self.uncertainties

    def _chi2(self, params: np.ndarray) -> float:
        """Compute χ²."""
        model_values = self.model(self.x_data, *params)
        return float(np.sum(((self.y_data - model_values) / self.uncertainties) ** 2))

    def fit(
        self,
        x0: Sequence[float],
        bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
    ) -> FitResult:
        """Perform χ² minimization.

        Parameters
        ----------
        x0 : sequence of float
            Initial parameter values.
        bounds : (lower, upper), optional
            Parameter bounds.

        Returns
        -------
        FitResult
        """
        x0_array = np.asarray(x0, dtype=np.float64)
        self._n_params = len(x0_array)

        ls_kwargs: dict[str, Any] = {
            "fun": self._residuals,
            "x0": x0_array,
            "method": self.config.method,
            "max_nfev": self.config.max_nfev,
            "x_scale": self.config.x_scale,
            "loss": self.config.loss,
            "f_scale": self.config.f_scale,
            "ftol": self.config.ftol,
            "xtol": self.config.xtol,
            "gtol": self.config.gtol,
        }
        if bounds is not None:
            ls_kwargs["bounds"] = (
                np.asarray(bounds[0], dtype=np.float64),
                np.asarray(bounds[1], dtype=np.float64),
            )

        result = least_squares(**ls_kwargs)

        # Compute covariance
        n = len(result.x)
        cov = np.zeros((n, n))
        if result.jac is not None:
            jtj = result.jac.T @ result.jac
            if np.linalg.cond(jtj) < 1e12:
                cov = np.linalg.inv(jtj)
            else:
                cov = np.linalg.pinv(jtj)
            chi2 = self._chi2(result.x)
            ndof = max(len(self.x_data) - n, 1)
            cov = cov * chi2 / ndof

        model_values = self.model(self.x_data, *result.x)

        fit_result = FitResult(
            parameters=result.x,
            covariance=cov,
            chi2=self._chi2(result.x),
            n_points=len(self.x_data),
            n_free=self._n_params,
            model_values=model_values,
            residuals=self.y_data - model_values,
            success=result.success,
            message=result.message,
            _residual_fn=self._residuals,
            _x_data=self.x_data,
            _y_data=self.y_data,
            _uncertainties=self.uncertainties,
        )

        self._last_result = fit_result
        return fit_result

    @property
    def parameters(self) -> np.ndarray:
        """Current best-fit parameters."""
        if self._last_result is None:
            raise RuntimeError("No fit result. Call fit() first.")
        return self._last_result.parameters

    @property
    def covariance(self) -> np.ndarray:
        """Current covariance matrix."""
        if self._last_result is None:
            raise RuntimeError("No fit result. Call fit() first.")
        return self._last_result.covariance
