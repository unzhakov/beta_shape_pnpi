"""
C(W) shape factor and g_V/g_A extraction from fitted spectra.

Provides methods to:
- Extract C(W) from fitted spectrum vs measured spectrum
- Perform Kurie plot analysis for endpoint determination
- Fit C(W) parametrization to extract g_V, g_A

The beta decay spectrum is:
    dΓ/dE = C₀ · p·E·(W₀-E)² · F(Z,E) · C(W) · [1 + δ(E)]

For allowed transitions with C(W) = 1, the Kurie plot
K(E) = sqrt(dΓ/dE / (p·E·F(Z,E)·[1+δ(E)]))
is linear in E with slope -(W₀-E₀)/E₀.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Optional, Tuple

import numpy as np

if TYPE_CHECKING:
    from beta_spectrum.spectrum import SpectrumConfig

from beta_spectrum.utils import T_to_W, momentum


@dataclass
class CWExtractionResult:
    """Container for C(W) extraction results."""

    energies_W: np.ndarray
    """Total energy grid in m_e units."""

    energies_keV: np.ndarray
    """Kinetic energy grid in keV."""

    cw_values: np.ndarray
    """Extracted C(W) values."""

    cw_errors: np.ndarray
    """Uncertainties on C(W) values."""

    flux: np.ndarray
    """Measured flux (counts per energy bin)."""

    flux_errors: np.ndarray
    """Uncertainties on flux."""

    endpoint_W: float
    """Endpoint energy W₀ in m_e units."""

    endpoint_keV: float
    """Endpoint kinetic energy in keV."""

    fit_result: Optional[Any] = None  # FitResult from SpectrumFitter

    def kurie_plot(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute Kurie plot values.

        Returns
        -------
        (energies_keV, kurie_values)
        """
        p = momentum(self.energies_W)
        E = self.energies_W

        # C(W) = (dΓ/dE) / (C₀ · p·E·(W₀-E)² · F · corrections)
        # For Kurie plot: K(E) = sqrt(C(W) / (p·E·(W₀-E)²))
        # Simplified: K(E) = sqrt(counts / (p·E·(W₀-E)²))
        phase_space = p * E * (self.endpoint_W - E) ** 2
        mask = phase_space > 0

        kurie = np.zeros_like(self.energies_keV)
        kurie[mask] = np.sqrt(
            self.flux[mask] / (phase_space[mask] * self.energies_keV[mask] ** 2)
        )

        return self.energies_keV[mask], kurie[mask]


@dataclass
class GVAExtractionResult:
    """Container for g_V and g_A extraction results."""

    gv_eff: float
    """Effective vector coupling constant."""

    gv_error: float
    """1σ uncertainty on g_V."""

    ga_eff: float
    """Effective axial coupling constant."""

    ga_error: float
    """1σ uncertainty on g_A."""

    correlation: float
    """Correlation coefficient between g_V and g_A."""

    chi2_per_dof: float
    """Reduced χ²."""

    def summary(self) -> str:
        lines = [
            "=" * 60,
            "g_V / g_A Extraction",
            "=" * 60,
            f"g_V^eff = {self.gv_eff:.6f} ± {self.gv_error:.6f}",
            f"g_A^eff = {self.ga_eff:.6f} ± {self.ga_error:.6f}",
            f"Correlation(g_V, g_A) = {self.correlation:.3f}",
            f"χ²/ndof = {self.chi2_per_dof:.3f}",
            "=" * 60,
        ]
        return "\n".join(lines)


class CWExtractor:
    """Extract C(W) shape factor from fitted spectrum.

    Parameters
    ----------
    config : SpectrumConfig
        Configuration used for the theoretical model.
    """

    def __init__(self, config: "SpectrumConfig"):
        self.config = config

    def extract(
        self,
        measured_counts: np.ndarray,
        measured_errors: np.ndarray,
        energies_keV: np.ndarray,
        model_counts: np.ndarray,
        endpoint_keV: float,
    ) -> CWExtractionResult:
        """Extract C(W) from measured spectrum and model prediction.

        Parameters
        ----------
        measured_counts : np.ndarray
            Measured counts per channel.
        measured_errors : np.ndarray
            Measurement uncertainties.
        energies_keV : np.ndarray
            Energy grid in keV.
        model_counts : np.ndarray
            Theoretical model counts (from fit).
        endpoint_keV : float
            Endpoint energy in keV.

        Returns
        -------
        CWExtractionResult
        """
        W = T_to_W(energies_keV)
        W0 = T_to_W(endpoint_keV)

        # C(W) ∝ measured / model (ratio gives shape factor)
        # For allowed transitions, model includes phase space + Fermi + corrections
        cw_values = np.where(model_counts > 0, measured_counts / model_counts, 1.0)
        cw_errors = np.where(
            model_counts > 0,
            measured_errors / model_counts,
            0.0,
        )

        # Clip C(W) to reasonable range
        cw_values = np.clip(cw_values, 0.01, 10.0)
        # Ensure errors are positive
        cw_errors = np.maximum(cw_errors, 1e-10)

        return CWExtractionResult(
            energies_W=W,
            energies_keV=energies_keV,
            cw_values=cw_values,
            cw_errors=cw_errors,
            flux=measured_counts,
            flux_errors=measured_errors,
            endpoint_W=W0,
            endpoint_keV=endpoint_keV,
        )

    def fit_parametrization(
        self,
        cw_result: "CWExtractionResult",
        parametrization: str = "constant",
    ) -> np.ndarray:
        """Fit C(W) parametrization to extracted values.

        Parameters
        ----------
        cw_result : CWExtractionResult
            Extracted C(W) values.
        parametrization : str
            Parametrization type: 'constant', 'linear', 'quadratic'.

        Returns
        -------
        np.ndarray
            Fitted parameters [a0, a1, ...].

        Raises
        ------
        ValueError
            If parametrization is not recognized.
        """
        valid = cw_result.energies_W > 0
        W = cw_result.energies_W[valid] - 1.0
        cw = cw_result.cw_values[valid]

        if len(W) < 3:
            raise ValueError("Not enough points for fit")

        if parametrization == "constant":
            coeffs, cov = np.polyfit(W, cw, 0, cov=True)
            return np.array([float(coeffs)])
        elif parametrization == "linear":
            coeffs, cov = np.polyfit(W, cw, 1, cov=True)
            return coeffs
        elif parametrization == "quadratic":
            coeffs, cov = np.polyfit(W, cw, 2, cov=True)
            return coeffs
        else:
            raise ValueError(f"Unknown parametrization: {parametrization}")

    def fit_gV_gA(
        self,
        cw_result: "CWExtractionResult",
        M_F: float = 1.0,
        M_GT: float = 0.0,
    ) -> "GVAExtractionResult":
        """Extract effective g_V and g_A from C(W) values.

        Parameters
        ----------
        cw_result : CWExtractionResult
            Extracted C(W) values.
        M_F : float
            Fermi matrix element.
        M_GT : float
            Gamow-Teller matrix element.

        Returns
        -------
        GVAExtractionResult
        """
        return GVAExtractor.extract_from_curve(
            cw_values=cw_result.cw_values,
            energies_W=cw_result.energies_W,
            endpoint_W=cw_result.endpoint_W,
        )

    def kurie_analysis(
        self,
        endpoint_keV: float,
    ) -> Tuple["CWExtractionResult", Optional[Any]]:
        """Perform Kurie plot analysis to extract C(W).

        Placeholder method: in a full implementation, this would
        use the fitted model to compute C(W) from experimental data.

        Returns
        -------
        (CWExtractionResult, fit_result)
            Extracted C(W) and optional fit result.
        """
        # For now, return a dummy result — the real implementation
        # would use the fitted spectrum model
        energies_keV = np.linspace(5, endpoint_keV * 0.98, 50)
        W = T_to_W(energies_keV / 1000.0)
        cw_values = np.ones_like(W)
        cw_errors = np.ones_like(W) * 0.05
        flux = np.ones_like(W) * 10.0
        flux_errors = np.ones_like(W) * 3.0
        W0 = T_to_W(endpoint_keV / 1000.0)

        result = CWExtractionResult(
            energies_W=W,
            energies_keV=energies_keV,
            cw_values=cw_values,
            cw_errors=cw_errors,
            flux=flux,
            flux_errors=flux_errors,
            endpoint_W=W0,
            endpoint_keV=endpoint_keV,
        )
        return result, None


class GVAExtractor:
    """Extract effective g_V and g_A from C(W) shape factor.

    For allowed transitions:
        C(W) ≈ 1 + (g_A/g_V)² · (corrections)

    The ratio g_A^eff/g_V^eff can be extracted from the shape of C(W).
    """

    @staticmethod
    def extract_from_curve(
        cw_values: np.ndarray,
        energies_W: np.ndarray,
        endpoint_W: float,
        transition_type: str = "F",
    ) -> GVAExtractionResult:
        """Extract g_V and g_A from C(W) values.

        Parameters
        ----------
        cw_values : np.ndarray
            Extracted C(W) values.
        energies_W : np.ndarray
            Total energy grid.
        endpoint_W : float
            Endpoint energy W₀.
        transition_type : str
            Transition type (F, GT, etc.).

        Returns
        -------
        GVAExtractionResult
        """
        # For allowed transitions, C(W) ≈ 1 + b·(W - W₀) where b depends on g_A/g_V
        # Simple linear fit of C(W) vs energy near endpoint
        mask = (energies_W > 0) & (energies_W < endpoint_W) & (cw_values > 0)

        if np.sum(mask) < 3:
            return GVAExtractionResult(
                gv_eff=1.0,
                gv_error=0.0,
                ga_eff=0.0,
                ga_error=0.0,
                correlation=0.0,
                chi2_per_dof=0.0,
            )

        # Fit C(W) = a + b·(W - W₀)
        x = energies_W[mask] - endpoint_W
        y = cw_values[mask]

        # Linear fit
        coeffs, cov = np.polyfit(x, y, 1, cov=True)
        a, b = coeffs

        # For allowed transitions: g_V ≈ 1 (CVC), g_A extracted from slope
        gv_eff = 1.0
        gv_error = 0.0
        ga_eff = b * endpoint_W  # Simplified relation
        ga_error = np.sqrt(cov[1, 1]) * endpoint_W

        correlation = (
            cov[0, 1] / (np.sqrt(cov[0, 0]) * np.sqrt(cov[1, 1]))
            if cov[0, 0] * cov[1, 1] > 0
            else 0.0
        )

        # χ² estimate
        y_fit = a + b * x
        chi2 = float(np.sum(((y - y_fit) / np.maximum(np.abs(y) * 0.01, 1e-10)) ** 2))
        ndof = max(len(x) - 2, 1)

        return GVAExtractionResult(
            gv_eff=gv_eff,
            gv_error=gv_error,
            ga_eff=ga_eff,
            ga_error=ga_error,
            correlation=float(correlation),
            chi2_per_dof=chi2 / ndof if ndof > 0 else chi2,
        )
