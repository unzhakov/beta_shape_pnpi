# beta_spectrum/fitter.py

"""
χ² minimization framework for beta spectrum fitting.

Provides infrastructure for fitting theoretical spectra to experimental data,
including parameter optimization, covariance estimation, and goodness-of-fit
metrics.

Usage:
    fitter = CurveFitter(model, data, uncertainties)
    result = fitter.fit(x0=[...], bounds=[...])
    print(f"Best fit: {result.parameters}")
    print(f"Covariance: {result.covariance}")
    print(f"χ²/ndf: {result.chi2_per_dof:.3f}")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.optimize import least_squares


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
    """Data minus model at best-fit parameters."""

    success: bool
    """Whether optimization converged."""

    message: str
    """Optimization success/failure message."""

    _model: Optional[Callable] = field(default=None, repr=False)
    """Model function (for profile likelihood)."""

    _x_data: Optional[np.ndarray] = field(default=None, repr=False)
    """Independent variable data (for profile likelihood)."""

    _y_data: Optional[np.ndarray] = field(default=None, repr=False)
    """Dependent variable data (for profile likelihood)."""

    _uncertainties: Optional[np.ndarray] = field(default=None, repr=False)
    """Measurement uncertainties (for profile likelihood)."""

    @property
    def chi2_per_dof(self) -> float:
        """χ² per degree of freedom."""
        return self.chi2 / self.n_free

    @property
    def p_value(self) -> float:
        """
        p-value for the fit (survival function of χ² distribution).

        High p-value (>0.05) suggests good fit; low p-value (<0.01) suggests
        poor fit or underestimated uncertainties.
        """
        from scipy.stats import chi2

        return 1.0 - chi2.cdf(self.chi2, self.n_points - self.n_free)

    @property
    def parameters_with_errors(self) -> Dict[str, Tuple[float, float]]:
        """
        Return parameter names with values and 1σ uncertainties.

        Returns dict mapping parameter names to (value, error) tuples.
        """
        errors = np.sqrt(np.diag(self.covariance))
        return {
            f"p{i}": (float(self.parameters[i]), float(errors[i]))
            for i in range(len(self.parameters))
        }

    @property
    def correlation_matrix(self) -> np.ndarray:
        """Parameter correlation matrix."""
        diag = np.sqrt(np.diag(self.covariance))
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
            f"χ²/ndof = {self.chi2_per_dof:.3f} " f"(ndof = {self.n_free})",
            f"p-value = {self.p_value:.4f}",
            "-" * 60,
            "Parameters:",
        ]

        errors = np.sqrt(np.diag(self.covariance))
        for i in range(len(self.parameters)):
            name = param_names[i] if param_names and i < len(param_names) else f"p{i}"
            lines.append(
                f"  {name:15s} = {self.parameters[i]:12.6f} " f"± {errors[i]:10.6f}"
            )

        lines.append("-" * 60)
        lines.append("Correlation matrix (off-diagonal):")
        corr = self.correlation_matrix
        n = len(self.parameters)
        header = "          " + "".join(f"{f'p{i}':>10s}" for i in range(n))
        lines.append(header)
        for i in range(n):
            row = f"  p{i!s:<8s}" + "".join(
                f"{corr[i, j]:10.3f}" for j in range(n) if j != i
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
        """
        Compute profile likelihood for a single parameter.

        Parameters
        ----------
        param_index : int
            Index of parameter to profile over.
        n_points : int
            Number of points to evaluate.
        param_range : (low, high), optional
            Range to scan. Uses ±3σ from fit result if None.

        Returns
        -------
        param_values : np.ndarray
            Parameter values scanned.
        chi2_values : np.ndarray
            χ² values at each point.
        normalized : np.ndarray
            χ² normalized to minimum (Δχ² = 0 at best fit).
        """
        if self._model is None or self._x_data is None:
            raise RuntimeError(
                "Cannot compute profile likelihood: model or data not stored. "
                "Use CurveFitter.profile_likelihood() instead."
            )

        if param_range is None:
            sigma = np.sqrt(abs(self.covariance[param_index, param_index]))
            center = self.parameters[param_index]
            param_range = (center - 3 * sigma, center + sigma * 3)

        param_values = np.linspace(param_range[0], param_range[1], n_points)
        chi2_values = np.zeros(n_points)

        unc = (
            self._uncertainties
            if self._uncertainties is not None
            else np.ones_like(self._x_data)
        )

        for i, pval in enumerate(param_values):
            params = list(self.parameters)
            params[param_index] = pval
            y_model = self._model(self._x_data, *params)
            residuals = (self._y_data - y_model) / unc
            chi2_values[i] = np.sum(residuals**2)

        min_chi2 = np.min(chi2_values)
        normalized = chi2_values - min_chi2

        return param_values, chi2_values, normalized

    def confidence_interval(
        self,
        param_index: int,
        confidence: float = 0.95,
        n_points: int = 100,
    ) -> Tuple[float, float]:
        """
        Compute confidence interval for a parameter using Δχ² method.

        Parameters
        ----------
        param_index : int
            Index of parameter.
        confidence : float
            Confidence level (e.g., 0.95 for 95% CI).
        n_points : int
            Number of points for profile likelihood scan.

        Returns
        -------
        (lower, upper) : float tuple
            Confidence interval bounds.
        """
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


@dataclass
class FitConfig:
    """Configuration for curve fitting."""

    method: str = "trf"
    """Optimization method passed to scipy.optimize.least_squares.
    Options: 'trf', 'dogbox', 'lm'."""

    max_nfev: int = 1000
    """Maximum number of function evaluations."""

    x_scale: Union[str, float] = "jac"
    """Parameter scaling: 'jac', 'diag', or float."""

    loss: str = "linear"
    """Loss function: 'linear' (χ²), 'soft_l1', 'huber', 'cauchy'."""

    f_scale: float = 1.0
    """Parameter for robust loss functions."""

    ftol: float = 1e-8
    """Termination tolerance for function value."""

    xtol: float = 1e-8
    """Termination tolerance for parameter change."""

    gtol: float = 1e-8
    """Termination tolerance for gradient."""


class CurveFitter:
    """
    χ² minimization framework for fitting models to data.

    The model function should have signature:
        model(x, *params) -> y

    where x is the independent variable and params are the fit parameters.

    χ² is computed as:
        χ² = Σ [(y_data[i] - model(x[i], params)) / σ[i]]²

    Parameters are optimized using scipy.optimize.least_squares with
    configurable loss function for robust fitting.
    """

    def __init__(
        self,
        model: Callable,
        x_data: np.ndarray,
        y_data: np.ndarray,
        uncertainties: Optional[np.ndarray] = None,
        config: Optional[FitConfig] = None,
    ):
        """
        Initialize fitter.

        Parameters
        ----------
        model : callable
            Model function: model(x, *params) -> y
        x_data : np.ndarray
            Independent variable data points.
        y_data : np.ndarray
            Dependent variable data points.
        uncertainties : np.ndarray, optional
            Measurement uncertainties (σ). If None, assumed to be 1.0.
        config : FitConfig, optional
            Fitting configuration. Uses defaults if None.
        """
        self.x_data = np.asarray(x_data, dtype=np.float64)
        self.y_data = np.asarray(y_data, dtype=np.float64)
        self.model = model

        if uncertainties is not None:
            self.uncertainties = np.asarray(uncertainties, dtype=np.float64)
        else:
            self.uncertainties = np.ones_like(self.y_data)

        # Filter out points with zero or negative uncertainties
        valid = self.uncertainties > 0
        self.x_data = self.x_data[valid]
        self.y_data = self.y_data[valid]
        self.uncertainties = self.uncertainties[valid]

        if config is not None:
            self.config = config
        else:
            self.config = FitConfig()

        self._n_params: Optional[int] = None

    def _residuals(self, params: np.ndarray) -> np.ndarray:
        """Compute weighted residuals for least_squares."""
        y_model = self.model(self.x_data, *params)
        return (self.y_data - y_model) / self.uncertainties

    def _chi2(self, params: np.ndarray) -> float:
        """Compute χ² for given parameters."""
        residuals = self._residuals(params)
        return float(np.sum(residuals**2))

    def fit(
        self,
        x0: Sequence[float],
        bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
        param_names: Optional[List[str]] = None,
    ) -> FitResult:
        """
        Perform χ² minimization.

        Parameters
        ----------
        x0 : sequence of float
            Initial parameter values.
        bounds : (lower, upper), optional
            Parameter bounds as (lower_bounds, upper_bounds).
            Use -np.inf / np.inf for unbounded parameters.
        param_names : list of str, optional
            Names for parameters (for summary output).

        Returns
        -------
        FitResult
            Fit results including parameters, covariance, χ², etc.
        """
        x0 = np.asarray(x0, dtype=np.float64)
        self._n_params = len(x0)

        if bounds is not None:
            lb = np.asarray(bounds[0], dtype=np.float64)
            ub = np.asarray(bounds[1], dtype=np.float64)
            bounds = (lb, ub)

        ls_kwargs = dict(
            fun=self._residuals,
            x0=x0,
            method=self.config.method,
            max_nfev=self.config.max_nfev,
            x_scale=self.config.x_scale,
            loss=self.config.loss,
            f_scale=self.config.f_scale,
            ftol=self.config.ftol,
            xtol=self.config.xtol,
            gtol=self.config.gtol,
        )
        if bounds is not None:
            ls_kwargs["bounds"] = bounds

        result = least_squares(**ls_kwargs)

        # Compute covariance matrix from Jacobian
        # Cov = (J^T J)^(-1) * χ² / (n - p)
        cov = self._compute_covariance(result)

        model_values = self.model(self.x_data, *result.x)
        residuals = self.y_data - model_values

        fit_result = FitResult(
            parameters=result.x,
            covariance=cov,
            chi2=self._chi2(result.x),
            n_points=len(self.x_data),
            n_free=self._n_params,
            model_values=model_values,
            residuals=residuals,
            success=result.success,
            message=result.message,
            _model=self.model,
            _x_data=self.x_data,
            _y_data=self.y_data,
            _uncertainties=self.uncertainties,
        )

        return fit_result

    def _compute_covariance(self, result) -> np.ndarray:
        """Estimate parameter covariance matrix from Jacobian."""
        n = len(result.x)
        cov = np.zeros((n, n))

        if result.jac is not None:
            # J^T J approximation
            jtj = result.jac.T @ result.jac
            if np.linalg.cond(jtj) < 1e12:
                cov = np.linalg.inv(jtj)
                # Scale by χ²/ndof
                chi2 = self._chi2(result.x)
                ndof = max(len(self.x_data) - n, 1)
                cov *= chi2 / ndof
            else:
                # Singular matrix — use pseudo-inverse
                cov = np.linalg.pinv(jtj)
                chi2 = self._chi2(result.x)
                ndof = max(len(self.x_data) - n, 1)
                cov *= chi2 / ndof

        return cov

    def residuals_plot(
        self,
        save_path: Optional[str] = None,
    ):
        """
        Plot data, model, and residuals.

        Parameters
        ----------
        save_path : str, optional
            Path to save figure.
        """
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(8, 8),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True,
        )

        # Data + model
        ax1.errorbar(
            self.x_data,
            self.y_data,
            yerr=self.uncertainties,
            fmt="ko",
            markersize=4,
            capsize=2,
            label="Data",
            alpha=0.7,
        )
        ax1.plot(
            self.x_data,
            self.model_values,
            "r-",
            lw=1.5,
            label="Model",
        )
        ax1.set_ylabel("Counts")
        ax1.set_title("Data vs. Model")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Residuals
        residuals = self.residuals / self.uncertainties
        ax2.plot(self.x_data, residuals, "bo", markersize=4, alpha=0.7)
        ax2.axhline(y=0, color="k", linestyle="--", lw=0.5)
        ax2.set_xlabel("Energy [MeV]")
        ax2.set_ylabel("Residuals (σ)")
        ax2.set_title(f"Residuals  χ²/ndof = {self.chi2_per_dof:.3f}")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")

        plt.show()

    @property
    def parameters(self) -> np.ndarray:
        """Current best-fit parameters (set by fit())."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError("No fit result available. Call fit() first.")
        return self._last_result.parameters

    @property
    def covariance(self) -> np.ndarray:
        """Current covariance matrix (set by fit())."""
        if not hasattr(self, "_last_result"):
            raise RuntimeError("No fit result available. Call fit() first.")
        return self._last_result.covariance

    def fit_and_store(
        self,
        x0: Sequence[float],
        bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
        param_names: Optional[List[str]] = None,
    ) -> FitResult:
        """
        Fit and store result for property access.

        Parameters
        ----------
        x0 : sequence of float
            Initial parameter values.
        bounds : (lower, upper), optional
            Parameter bounds.
        param_names : list of str, optional
            Parameter names for summary.

        Returns
        -------
        FitResult
        """
        result = self.fit(x0, bounds, param_names)
        self._last_result = result
        return result
