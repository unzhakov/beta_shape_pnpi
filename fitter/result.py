"""
Analysis result container for beta spectrum fitting.

Wraps fit results with physics-specific metadata (nuclide, extracted
C(W) values, g_V, g_A, etc.).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

import numpy as np

if TYPE_CHECKING:
    from beta_spectrum.spectrum import SpectrumConfig


@dataclass
class FitSummary:
    """Compact fit summary for reporting."""

    nuclide: str
    """Nuclide identifier (e.g. 'Tc99')."""

    converged: bool
    """Whether the fit converged."""

    chi2_per_dof: float
    """Reduced χ²."""

    p_value: float
    """Fit p-value."""

    parameters: Dict[str, Tuple[float, float]]
    """Parameter names → (value, error)."""

    def __str__(self) -> str:
        lines = [
            f"Fit Summary: {self.nuclide}",
            f"  Converged: {self.converged}",
            f"  χ²/ndof: {self.chi2_per_dof:.3f}",
            f"  p-value: {self.p_value:.4f}",
            "  Parameters:",
        ]
        for name, (val, err) in self.parameters.items():
            lines.append(f"    {name:20s} = {val:12.6f} ± {err:10.6f}")
        return "\n".join(lines)


@dataclass
class AnalysisResult:
    """Complete analysis result from spectrum fitting.

    Parameters
    ----------
    fit_result : FitResult
        Low-level fit result from SpectrumFitter.
    nuclide : str, optional
        Nuclide identifier.
    config : SpectrumConfig, optional
        Original spectrum configuration.
    cw_values : np.ndarray, optional
        Extracted C(W) values at each energy bin.
    cw_errors : np.ndarray, optional
        Uncertainties on C(W).
    gv_eff : float, optional
        Effective vector coupling constant.
    ga_eff : float, optional
        Effective axial coupling constant.
    endpoint_keV : float, optional
        Extracted endpoint energy.
    metadata : dict, optional
        Additional analysis metadata.
    """

    fit_result: Any  # FitResult from SpectrumFitter
    nuclide: str = ""
    config: Optional["SpectrumConfig"] = None
    cw_values: Optional[np.ndarray] = None
    cw_errors: Optional[np.ndarray] = None
    gv_eff: Optional[float] = None
    ga_eff: Optional[float] = None
    endpoint_keV: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def summary(self) -> FitSummary:
        """Get a FitSummary from this result."""
        fr = self.fit_result
        param_names = self.metadata.get("param_names", None)
        params = {}
        if hasattr(fr, "parameters_with_errors"):
            params = fr.parameters_with_errors
        elif hasattr(fr, "parameters"):
            errors = np.sqrt(np.maximum(np.diag(fr.covariance), 0.0))
            for i in range(len(fr.parameters)):
                name = (
                    param_names[i] if param_names and i < len(param_names) else f"p{i}"
                )
                params[name] = (float(fr.parameters[i]), float(errors[i]))

        return FitSummary(
            nuclide=self.nuclide,
            converged=fr.success if hasattr(fr, "success") else False,
            chi2_per_dof=fr.chi2_per_dof if hasattr(fr, "chi2_per_dof") else 0.0,
            p_value=fr.p_value if hasattr(fr, "p_value") else 0.0,
            parameters=params,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "nuclide": self.nuclide,
            "converged": (
                self.fit_result.success
                if hasattr(self.fit_result, "success")
                else False
            ),
            "chi2_per_dof": (
                self.fit_result.chi2_per_dof
                if hasattr(self.fit_result, "chi2_per_dof")
                else 0.0
            ),
            "p_value": (
                self.fit_result.p_value if hasattr(self.fit_result, "p_value") else 0.0
            ),
            "gv_eff": self.gv_eff,
            "ga_eff": self.ga_eff,
            "endpoint_keV": self.endpoint_keV,
            "metadata": self.metadata,
        }
