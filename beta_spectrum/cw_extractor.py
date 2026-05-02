# mypy: disable-error-code = no-any-return, arg-type, assignment, call-arg, return-value, no-untyped-def, misc
# beta_spectrum/cw_extractor.py
# ruff: noqa: ANN201, ANN001, ANN204
"""
C(W) shape factor extraction from beta spectra.

Provides methods to:
- Extract C(W) from measured beta spectra
- Perform Kurie plot analysis for endpoint determination
- Fit C(W) parametrization to extract g_V, g_A
- Compare extracted C(W) with theoretical predictions

The beta decay spectrum is:
    dΓ/dE = C₀ · p·E·(W₀-E)² · F(Z,E) · C(W) · [1 + δ(E)]

where:
    p, E  = momentum, total energy of electron
    W₀    = endpoint total energy in m_e units
    F(Z,E) = Fermi function
    C(W)  = shape factor (what we extract)
    δ(E)  = radiative + screening + other corrections

For allowed transitions with C(W) = 1, the Kurie plot
K(E) = sqrt(dΓ/dE / (p·E·F(Z,E)·[1+δ(E)]))
is linear in E with slope -(W₀-E₀)/E₀.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from beta_spectrum.components.detector_response import DetectorResponse

import numpy as np

from .fitter import CurveFitter, FitResult
from .spectrum import SpectrumConfig
from .utils import T_to_W, W_to_T, momentum


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

    fit_result: Optional[FitResult] = None
    """Fit result from Kurie plot endpoint determination."""

    def kurie_plot(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute Kurie plot values.

        For allowed transitions (C(W)=1), K(E) should be linear.
        Deviations from linearity indicate shape factor effects.
        """
        # K(E) = sqrt(flux / (p * E * F * phase_space * corrections))
        # We compute this from the extracted C(W) values
        # K(E) ∝ (W0 - E) / sqrt(C(W))
        kurie = np.sqrt(self.flux) / (
            momentum(self.energies_W)
            * np.sqrt(self.energies_W)
            * (self.endpoint_W - self.energies_W)
        )
        kurie = kurie / np.max(kurie)  # Normalize
        return self.energies_keV, kurie


@dataclass
class GVAExtractionResult:
    """Container for g_V, g_A extraction results."""

    g_V: float
    """Extracted vector coupling constant."""

    g_V_error: float
    """Uncertainty on g_V."""

    g_A: float
    """Extracted axial-vector coupling constant."""

    g_A_error: float
    """Uncertainty on g_A."""

    fit_result: FitResult
    """Fit result with full covariance matrix."""

    cw_data: np.ndarray
    """C(W) data used in fit."""

    cw_errors: np.ndarray
    """C(W) uncertainties."""

    energy_W: np.ndarray
    """Energy grid for C(W) data."""

    def summary(self) -> str:
        """Human-readable summary of g_V, g_A extraction."""
        lines = [
            "=" * 60,
            "g_V and g_A Extraction Results",
            "=" * 60,
            f"g_V  = {self.g_V:.4f} ± {self.g_V_error:.4f}",
            f"g_A  = {self.g_A:.4f} ± {self.g_A_error:.4f}",
            f"χ²/ndf = {self.fit_result.chi2_per_dof:.3f}",
            f"p-value = {self.fit_result.p_value:.4f}",
            "-" * 60,
            "Comparison with standard values:",
            f"  g_V / g_A ratio = {self.g_V / self.g_A:.4f} "
            f"(± {self._ratio_error():.4f})",
            "  Standard model: g_V/g_A ≈ 1.0 for pure Fermi",
            "=" * 60,
        ]
        return "\n".join(lines)

    def _ratio_error(self) -> float:
        """Error on g_V/g_A ratio using error propagation."""
        r = self.g_V / self.g_A
        cov = self.fit_result.covariance
        dg_V = np.sqrt(cov[0, 0])
        dg_A = np.sqrt(cov[1, 1])
        dcorr = cov[0, 1]
        return abs(r) * np.sqrt(
            (dg_V / self.g_V) ** 2
            + (dg_A / self.g_A) ** 2
            - 2 * dcorr / (self.g_V * self.g_A)
        )


class CWExtractor:
    """
    Extract C(W) shape factor from measured beta spectra.

    The extraction works by dividing the measured spectrum by all known
    theoretical factors:

        C(W) = (dΓ/dE)_measured / [C₀ · p·E·(W₀-E)² · F(Z,E) · (1+δ)]

    where C₀ is an overall normalization constant.

    For allowed transitions, C(W) should be approximately constant (=1).
    For forbidden transitions, C(W) has energy dependence that encodes
    nuclear structure information.
    """

    def __init__(
        self,
        config: SpectrumConfig,
        flux: np.ndarray,
        energies_keV: np.ndarray,
        flux_errors: Optional[np.ndarray] = None,
        detector_response: Optional["DetectorResponse"] = None,  # noqa: F821
    ):
        """
        Initialize extractor.

        Parameters
        ----------
        config : SpectrumConfig
            Decay configuration (Z, A, endpoint, corrections).
        flux : np.ndarray
            Measured spectrum (counts per energy bin).
        energies_keV : np.ndarray
            Kinetic energy grid in keV.
        flux_errors : np.ndarray, optional
            Uncertainties on flux. If None, assumed to be sqrt(flux).
        detector_response : DetectorResponse, optional
            Detector response for deconvolution (optional).
        """
        self.config = config
        self.flux = np.asarray(flux, dtype=np.float64)
        self.energies_keV = np.asarray(energies_keV, dtype=np.float64)
        self.energies_W = np.array([T_to_W(e / 1000.0) for e in self.energies_keV])
        self.detector_response = detector_response

        if flux_errors is not None:
            self.flux_errors = np.asarray(flux_errors, dtype=np.float64)
        else:
            # Poisson statistics
            self.flux_errors = np.sqrt(np.maximum(self.flux, 1.0))

    def _theoretical_factors(self, W: np.ndarray) -> Dict[str, np.ndarray]:
        """Compute all known theoretical factors."""
        factors = {}

        # Phase space: p * E
        p = momentum(W)
        factors["phase_space"] = p * W

        # Fermi function
        from .components.fermi import FermiFunction

        fermi = FermiFunction(self.config.Z_parent, self.config.A_number)
        factors["fermi"] = fermi(W)

        # Screening
        if self.config.use_screening:
            from .components.screening import ScreeningCorrection

            screening = ScreeningCorrection(
                FermiFunction(self.config.Z_parent, self.config.A_number)
            )
            factors["screening"] = screening(W)
        else:
            factors["screening"] = np.ones_like(W)

        # Exchange
        if self.config.use_exchange:
            from .components.exchange import ExchangeCorrection

            exchange = ExchangeCorrection(Z=self.config.Z_parent)
            factors["exchange"] = exchange(W)
        else:
            factors["exchange"] = np.ones_like(W)

        # Finite size
        if self.config.use_finite_size:
            from .components.finite_size import (
                FiniteSizeL0,
                ChargeDistributionU,
            )

            fs_l0 = FiniteSizeL0(Z=self.config.Z_daughter, A=self.config.A_number)
            fs_u = ChargeDistributionU(Z=self.config.Z_daughter, A=self.config.A_number)
            factors["finite_size_L0"] = fs_l0(W)
            factors["charge_distribution"] = fs_u(W)
        else:
            factors["finite_size_L0"] = np.ones_like(W)
            factors["charge_distribution"] = np.ones_like(W)

        # Radiative corrections
        if self.config.use_radiative:
            from .components.radiative import RadiativeCorrection

            W0 = T_to_W(self.config.endpoint_MeV)
            rad = RadiativeCorrection(
                W0=W0,
                Z=self.config.Z_parent,
                A=self.config.A_number,
                use_endpoint_resummation=True,
                delta_cut=1e-3,
            )
            factors["radiative"] = rad(W)
        else:
            factors["radiative"] = np.ones_like(W)

        return factors

    def extract_CW(
        self,
        endpoint_keV: Optional[float] = None,
        normalize_method: str = "endpoint",
    ) -> CWExtractionResult:
        """
        Extract C(W) from measured spectrum.

        Parameters
        ----------
        endpoint_keV : float, optional
            Endpoint energy in keV. If None, determined from Kurie plot.
        normalize_method : str
            Method for normalization: 'endpoint' or 'fit'.

        Returns
        -------
        CWExtractionResult
            Extracted C(W) values with uncertainties.
        """
        if endpoint_keV is not None:
            W0 = T_to_W(endpoint_keV / 1000.0)
        else:
            W0 = self._determine_endpoint()

        # Compute theoretical factors (excluding C(W))
        factors = self._theoretical_factors(self.energies_W)

        # Product of all known factors
        product = factors["phase_space"] * factors["fermi"]
        for key in [
            "screening",
            "exchange",
            "finite_size_L0",
            "charge_distribution",
            "radiative",
        ]:
            product *= factors[key]

        # Endpoint factor: (W0 - W)^2
        endpoint_factor = (W0 - self.energies_W) ** 2

        # Overall normalization
        if normalize_method == "endpoint":
            # Use flux at lowest energy where C(W) ≈ 1 for allowed decay
            low_E_mask = self.energies_W < W0 - 0.02
            if np.any(low_E_mask):
                C0 = np.median(
                    self.flux[low_E_mask]
                    / (product[low_E_mask] * endpoint_factor[low_E_mask])
                )
            else:
                C0 = 1.0
        else:
            # Fit normalization
            C0 = self._fit_normalization(W0, product, endpoint_factor)

        # Extract C(W)
        cw = self.flux / (C0 * product * endpoint_factor)

        # Propagate uncertainties
        cw_errors = np.maximum(np.abs(cw), 1e-10) * np.sqrt(
            (self.flux_errors / np.maximum(self.flux, 1e-30)) ** 2
            + 0.01**2  # ~1% systematic on theoretical factors
        )

        # Clip unphysical values
        cw = np.clip(cw, 0.01, 10.0)
        cw_errors = np.clip(cw_errors, 0.0, np.inf)

        return CWExtractionResult(
            energies_W=self.energies_W.copy(),
            energies_keV=self.energies_keV.copy(),
            cw_values=cw,
            cw_errors=cw_errors,
            flux=self.flux.copy(),
            flux_errors=self.flux_errors.copy(),
            endpoint_W=W0,
            endpoint_keV=W_to_T(W0) * 1000.0,
            fit_result=None,
        )

    def kurie_analysis(
        self,
        endpoint_keV: Optional[float] = None,
    ) -> Tuple[CWExtractionResult, FitResult]:
        """
        Perform Kurie plot analysis for endpoint determination.

        For allowed transitions, the Kurie plot K(E) = sqrt(Φ/(p·E·F))
        should be linear in E. The endpoint is where K(E) → 0.

        Parameters
        ----------
        endpoint_keV : float, optional
            If provided, use this endpoint. Otherwise fit for it.

        Returns
        -------
        cw_result : CWExtractionResult
            Extracted C(W) values.
        fit_result : FitResult
            Kurie plot fit result.
        """
        # First extract C(W) with a rough endpoint
        cw_result = self.extract_CW(
            endpoint_keV=endpoint_keV,
            normalize_method="fit",
        )

        if endpoint_keV is None:
            # Fit Kurie plot to determine endpoint
            fit_result = self._fit_kurie_endpoint(cw_result)
        else:
            # Dummy fit result
            fit_result = FitResult(
                parameters=np.array([cw_result.endpoint_keV]),
                covariance=np.eye(1) * 1.0,
                chi2=0.0,
                n_points=len(cw_result.energies_keV),
                n_free=1,
                model_values=np.zeros_like(cw_result.energies_keV),
                residuals=np.zeros_like(cw_result.energies_keV),
                success=True,
                message="Endpoint fixed, no fit needed",
            )

        cw_result.fit_result = fit_result

        return cw_result, fit_result

    def fit_CW_parametrization(
        self,
        cw_result: CWExtractionResult,
        parametrization: str = "linear",
        bounds: Optional[Tuple] = None,
    ) -> FitResult:
        """
        Fit C(W) to a parametrized form.

        Supported parametrizations:
        - "constant": C(W) = a₀ (allowed transition)
        - "linear": C(W) = a₀ + a₁·(W - 1) (first-order correction)
        - "quadratic": C(W) = a₀ + a₁·(W - 1) + a₂·(W - 1)²
        - "gV_gA": C(W) = g_V²·|M_F|² + g_A²·|M_GT|² + interference terms

        Parameters
        ----------
        cw_result : CWExtractionResult
            Extracted C(W) data.
        parametrization : str
            Parametrization form.
        bounds : tuple, optional
            Parameter bounds (lower, upper).

        Returns
        -------
        FitResult
            Fit results.
        """
        # Define parametrization models
        models = {
            "constant": self._model_constant,
            "linear": self._model_linear,
            "quadratic": self._model_quadratic,
        }

        if parametrization not in models:
            raise ValueError(
                f"Unknown parametrization: {parametrization}. "
                f"Use one of: {list(models.keys())}"
            )

        n_params = {
            "constant": 1,
            "linear": 2,
            "quadratic": 3,
        }[parametrization]

        # Initial guesses
        x0 = self._initial_guesses(parametrization, cw_result)

        # Filter out points near endpoint (large uncertainties)
        mask = cw_result.energies_keV < cw_result.endpoint_keV * 0.95
        energies = cw_result.energies_W[mask]
        cw = cw_result.cw_values[mask]
        cw_err = cw_result.cw_errors[mask]

        if len(energies) < n_params + 2:
            raise ValueError(
                f"Not enough data points for {parametrization} fit. "
                f"Need at least {n_params + 2}, got {len(energies)}"
            )

        fitter = CurveFitter(
            model=models[parametrization],
            x_data=energies,
            y_data=cw,
            uncertainties=cw_err,
        )

        return fitter.fit(x0, bounds, param_names=[f"a{i}" for i in range(n_params)])

    def fit_gV_gA(
        self,
        cw_result: CWExtractionResult,
        M_F: float = 1.0,
        M_GT: float = 1.0,
        interference: float = 0.0,
        bounds: Optional[Tuple] = None,
    ) -> GVAExtractionResult:
        """
        Extract g_V and g_A from C(W) data.

        For mixed transitions:
            C(W) = g_V²·|M_F|² + g_A²·|M_GT|² + b·g_V·g_A·interference

        where M_F and M_GT are Fermi and Gamow-Teller matrix elements.

        Parameters
        ----------
        cw_result : CWExtractionResult
            Extracted C(W) data.
        M_F : float
            Fermi matrix element (default: 1.0 for pure Fermi).
        M_GT : float
            Gamow-Teller matrix element (default: 1.0 for pure GT).
        interference : float
            Interference term coefficient (b coefficient).
        bounds : tuple, optional
            Parameter bounds (g_V_low, g_A_low, g_V_high, g_A_high).

        Returns
        -------
        GVAExtractionResult
            Extracted g_V, g_A with uncertainties.
        """
        # Filter points near endpoint
        mask = cw_result.energies_keV < cw_result.endpoint_keV * 0.95
        energies = cw_result.energies_W[mask]
        cw = cw_result.cw_values[mask]
        cw_err = cw_result.cw_errors[mask]

        if len(energies) < 4:
            raise ValueError(
                f"Not enough data points for g_V/g_A fit. "
                f"Need at least 4, got {len(energies)}"
            )

        # g_V, g_A model
        def gVgA_model(W, g_V, g_A):
            return g_V**2 * M_F**2 + g_A**2 * M_GT**2 + interference * g_V * g_A

        # Initial guesses: use Paulsen et al. values as starting point
        x0 = [0.376, 0.574]

        if bounds is None:
            bounds = (
                [0.0, 0.0],
                [2.0, 2.0],
            )

        fitter = CurveFitter(
            model=gVgA_model,
            x_data=energies,
            y_data=cw,
            uncertainties=cw_err,
        )

        fit_result = fitter.fit(
            x0,
            bounds,
            param_names=["g_V", "g_A"],
        )

        return GVAExtractionResult(
            g_V=float(fit_result.parameters[0]),
            g_V_error=float(np.sqrt(fit_result.covariance[0, 0])),
            g_A=float(fit_result.parameters[1]),
            g_A_error=float(np.sqrt(fit_result.covariance[1, 1])),
            fit_result=fit_result,
            cw_data=cw,
            cw_errors=cw_err,
            energy_W=energies,
        )

    def _determine_endpoint(self) -> float:
        """
        Determine endpoint energy from Kurie plot linearity.

        Scans endpoint values and finds the one that gives the most
        linear Kurie plot (highest R² for linear fit).
        """
        # Scan endpoint in keV
        E_min = self.energies_keV[-1] - 50  # At least 50 keV below max
        E_max = self.energies_keV[-1] + 10  # Allow slight overshoot

        best_R2 = -np.inf
        best_endpoint = E_max

        for E_trial in np.linspace(E_min, E_max, 30):
            W0 = T_to_W(E_trial / 1000.0)

            # Compute Kurie values
            mask = self.energies_W < W0 - 0.01
            if not np.any(mask):
                continue

            W = self.energies_W[mask]
            flux = self.flux[mask]

            # K(E) = sqrt(flux / (p*E*F))
            p = momentum(W)
            from .components.fermi import FermiFunction

            F = FermiFunction()(W, self.config.Z_parent)

            kurie = np.sqrt(flux / (p * W * F))

            # Linear fit: K(E) = a + b*E
            E_keV = self.energies_keV[mask]
            try:
                coeffs = np.polyfit(E_keV, kurie, 1)
                # R²
                kurie_fit = np.polyval(coeffs, E_keV)
                SS_res = np.sum((kurie - kurie_fit) ** 2)
                SS_tot = np.sum((kurie - np.mean(kurie)) ** 2)
                R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0

                if R2 > best_R2:
                    best_R2 = R2
                    best_endpoint = E_trial
            except np.linalg.LinAlgError:
                continue

        return T_to_W(best_endpoint / 1000.0)

    def _fit_kurie_endpoint(self, cw_result: CWExtractionResult) -> FitResult:
        """Fit Kurie plot to determine endpoint energy."""
        W0 = cw_result.endpoint_W

        # Kurie plot: K(E) = sqrt(flux / (p*E*F))
        mask = self.energies_W < W0 - 0.01
        W = self.energies_W[mask]
        flux = self.flux[mask]
        flux_err = self.flux_errors[mask]

        p = momentum(W)
        from .components.fermi import FermiFunction

        F = FermiFunction()(W, self.config.Z_parent)

        kurie = np.sqrt(flux / (p * W * F))
        E_keV = self.energies_keV[mask]

        # K(E) = a + b*(E - E0), where K(E0) = 0
        # So a = -b*(E0 - E0) = 0 at the true endpoint
        # Linear fit: K(E) = a + b*E
        def kurie_model(E, a, b):
            return a + b * E

        fitter = CurveFitter(
            model=kurie_model,
            x_data=E_keV,
            y_data=kurie,
            uncertainties=flux_err / (2 * kurie * p * W * F),
        )

        return fitter.fit(
            x0=[np.mean(kurie), -np.mean(kurie) / E_keV[-1]],
            param_names=["a", "b"],
        )

    def _fit_normalization(
        self, W0: float, product: np.ndarray, endpoint_factor: np.ndarray
    ) -> float:
        """Fit overall normalization C₀."""
        mask = self.energies_W < W0 - 0.02
        if not np.any(mask):
            return 1.0

        def norm_model(E, C0):
            return C0 * product[mask] * endpoint_factor[mask]

        fitter = CurveFitter(
            model=norm_model,
            x_data=self.energies_W[mask],
            y_data=self.flux[mask],
            uncertainties=self.flux_errors[mask],
        )

        result = fitter.fit(x0=[1.0])
        return float(result.parameters[0])

    def _model_constant(self, W: np.ndarray, a0: float) -> np.ndarray:
        """C(W) = a₀ (constant)."""
        return np.full_like(W, a0)

    def _model_linear(self, W: np.ndarray, a0: float, a1: float) -> np.ndarray:
        """C(W) = a₀ + a₁·(W - 1)."""
        return a0 + a1 * (W - 1.0)

    def _model_quadratic(
        self, W: np.ndarray, a0: float, a1: float, a2: float
    ) -> np.ndarray:
        """C(W) = a₀ + a₁·(W - 1) + a₂·(W - 1)²."""
        return a0 + a1 * (W - 1.0) + a2 * (W - 1.0) ** 2

    def _initial_guesses(
        self, parametrization: str, cw_result: CWExtractionResult
    ) -> List[float]:
        """Generate initial guesses for fit parameters."""
        # Average C(W) in mid-energy region
        mid_mask = (cw_result.energies_keV > cw_result.endpoint_keV * 0.2) & (
            cw_result.energies_keV < cw_result.endpoint_keV * 0.7
        )
        if np.any(mid_mask):
            avg_cw = float(np.mean(cw_result.cw_values[mid_mask]))
        else:
            avg_cw = 1.0

        if parametrization == "constant":
            return [avg_cw]
        elif parametrization == "linear":
            return [avg_cw, 0.0]
        elif parametrization == "quadratic":
            return [avg_cw, 0.0, 0.0]
        else:
            return [avg_cw] * 3


# mypy: ignore-errors
