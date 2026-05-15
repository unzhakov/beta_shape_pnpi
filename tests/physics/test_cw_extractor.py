# tests/test_cw_extractor.py — Tests for C(W) shape factor extraction

"""
Why test this?
--------------
The C(W) extractor is the core scientific tool for extracting the shape factor
from measured beta spectra. It must correctly handle:
 1. Energy unit conversions (keV ↔ m_e units)
 2. Theoretical factor computation (Fermi, screening, radiative, etc.)
 3. Normalization of extracted C(W)
 4. Kurie plot endpoint determination
 5. g_V, g_A extraction from C(W) data

We test with synthetic data where the true answer is known, so we can verify
the extraction recovers the correct values.
"""

import numpy as np
import pytest

from beta_spectrum import (
    BetaSpectrum,
    CWExtractor,
    CWExtractionResult,
    GVAExtractionResult,
    SpectrumConfig,
    T_to_W,
    W_to_T,
)


class TestCWExtractionBasic:
    """Basic C(W) extraction tests."""

    def _make_synthetic_data(
        self,
        endpoint_keV: float = 294.0,
        n_points: int = 100,
        seed: int = 42,
    ):
        """Create synthetic beta spectrum with C(W) = 1 (allowed decay)."""
        np.random.seed(seed)

        # Energy grid
        energies_keV = np.linspace(5, endpoint_keV * 0.98, n_points)

        # Build theoretical spectrum with C(W) = 1
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=endpoint_keV / 1000.0,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )
        spectrum = BetaSpectrum.from_config(config)
        W, E_MeV = spectrum.get_energy_grid(config)

        # Interpolate to our energy grid
        from scipy.interpolate import interp1d

        flux_interp = interp1d(
            E_MeV * 1000, spectrum(W), kind="linear", fill_value=0.0, bounds_error=False
        )
        flux_theory = flux_interp(energies_keV)

        # Add noise (Poisson-like)
        noise_std = np.sqrt(np.maximum(flux_theory, 1e-10))
        flux = flux_theory + np.random.normal(0, noise_std, len(flux_theory))
        flux = np.maximum(flux, 0.0)  # Flux cannot be negative
        flux_errors = np.sqrt(np.maximum(flux, 1e-10))

        return energies_keV, flux, flux_errors, config

    def test_CW_extraction_allowed_decay(self):
        """For allowed decay (C(W)=1), extracted C(W) should be ~1."""
        energies_keV, flux, flux_errors, config = self._make_synthetic_data()

        extractor = CWExtractor(config, flux, energies_keV, flux_errors)
        result = extractor.extract_CW(endpoint_keV=294.0)

        # C(W) should be close to 1 in mid-energy region
        mid_mask = (energies_keV > 50) & (energies_keV < 200)
        if np.any(mid_mask):
            avg_cw = np.mean(result.cw_values[mid_mask])
            assert (
                abs(avg_cw - 1.0) < 0.3
            ), f"Allowed decay C(W) should be ~1, got {avg_cw}"

    def test_CW_extraction_preserves_endpoint(self):
        """Extracted C(W) should report correct endpoint."""
        energies_keV, flux, flux_errors, config = self._make_synthetic_data()

        extractor = CWExtractor(config, flux, energies_keV, flux_errors)
        result = extractor.extract_CW(endpoint_keV=294.0)

        assert (
            abs(result.endpoint_keV - 294.0) < 1.0
        ), f"Endpoint should be ~294 keV, got {result.endpoint_keV}"

    def test_CW_extraction_errors_positive(self):
        """C(W) uncertainties should be positive."""
        energies_keV, flux, flux_errors, config = self._make_synthetic_data()

        extractor = CWExtractor(config, flux, energies_keV, flux_errors)
        result = extractor.extract_CW(endpoint_keV=294.0)

        assert np.all(result.cw_errors > 0), "All errors should be positive"

    def test_CW_extraction_clipped(self):
        """C(W) values should be clipped to [0.01, 10.0]."""
        energies_keV, flux, flux_errors, config = self._make_synthetic_data()

        extractor = CWExtractor(config, flux, energies_keV, flux_errors)
        result = extractor.extract_CW(endpoint_keV=294.0)

        assert np.all(result.cw_values >= 0.01), "C(W) should be >= 0.01"
        assert np.all(result.cw_values <= 10.0), "C(W) should be <= 10.0"


class TestCWParametrization:
    """Test C(W) parametrization fitting."""

    def _make_CW_data(
        self,
        endpoint_keV: float = 294.0,
        n_points: int = 60,
        seed: int = 42,
    ):
        """Create synthetic C(W) data with known parametrization."""
        np.random.seed(seed)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=endpoint_keV / 1000.0,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        # Generate C(W) values with known form
        W = np.linspace(T_to_W(0.05), T_to_W(endpoint_keV * 0.95), n_points)
        E_keV = np.array([W_to_T(w) for w in W])

        # Linear C(W): C(W) = 1.0 + 0.5*(W-1)
        true_a0, true_a1 = 1.0, 0.5
        cw_true = true_a0 + true_a1 * (W - 1.0)
        cw_noise = cw_true + np.random.normal(0, 0.05, len(W))

        # C(W) = flux / (C0 * phase_space * Fermi * corrections * (W0-W)^2)
        # So: flux = C(W) * C0 * phase_space * Fermi * corrections * (W0-W)^2

        p = np.sqrt(np.maximum(W**2 - 1.0, 0.0))
        from beta_spectrum.components.fermi import FermiFunction

        F = FermiFunction(Z=config.Z_parent, A=config.A_number)(W)

        C0 = 1.0  # Normalization
        endpoint_factor = (T_to_W(endpoint_keV / 1000.0) - W) ** 2

        flux_synthetic = cw_noise * C0 * p * W * F * endpoint_factor
        flux_synthetic = np.maximum(flux_synthetic, 0)
        flux_errors = np.sqrt(np.maximum(flux_synthetic, 1.0))

        return E_keV, flux_synthetic, flux_errors, config, true_a0, true_a1

    def test_constant_parametrization(self):
        """Fit C(W) = a₀ with known constant value."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw_true = np.full_like(W, 1.0)
        cw = cw_true + np.random.normal(0, 0.05, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.05
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_CW_parametrization(cw_result, parametrization="constant")

        assert (
            abs(result.parameters[0] - 1.0) < 0.1
        ), f"Constant C(W) should be ~1, got {result.parameters[0]}"

    def test_linear_parametrization(self):
        """Fit C(W) = a₀ + a₁·(W-1) with known coefficients."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        true_a0, true_a1 = 1.0, 0.5
        cw = true_a0 + true_a1 * (W - 1.0) + np.random.normal(0, 0.05, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.05
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_CW_parametrization(cw_result, parametrization="linear")

        assert (
            abs(result.parameters[0] - true_a0) < 0.15
        ), f"a₀: expected {true_a0}, got {result.parameters[0]}"
        assert (
            abs(result.parameters[1] - true_a1) < 0.5
        ), f"a₁: expected {true_a1}, got {result.parameters[1]}"

    def test_quadratic_parametrization(self):
        """Fit C(W) = a₀ + a₁·(W-1) + a₂·(W-1)² with known coefficients."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 60)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        true_a0, true_a1, true_a2 = 1.0, 0.3, -0.1
        cw = true_a0 + true_a1 * (W - 1.0) + true_a2 * (W - 1.0) ** 2
        cw += np.random.normal(0, 0.03, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.03
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.03,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.03,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_CW_parametrization(
            cw_result, parametrization="quadratic"
        )

        assert abs(result.parameters[0] - true_a0) < 0.15
        assert abs(result.parameters[1] - true_a1) < 0.5
        assert abs(result.parameters[2] - true_a2) < 0.5

    def test_invalid_parametrization(self):
        """Invalid parametrization should raise ValueError."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 30)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw = np.ones_like(W)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )

        extractor = CWExtractor(config, cw, energies_keV)
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        with pytest.raises(ValueError, match="Unknown parametrization"):
            extractor.fit_CW_parametrization(cw_result, parametrization="invalid")


class TestGVAExtraction:
    """Test g_V, g_A extraction from C(W) data."""

    def test_gV_gA_pure_fermi(self):
        """For pure Fermi transition, g_A ≈ 0."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])

        # Pure Fermi: C(W) = g_V² (constant)
        true_gV = 0.5
        cw = true_gV**2 + np.random.normal(0, 0.02, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.02
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.02,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.02,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_gV_gA(cw_result, M_F=1.0, M_GT=0.0)

        assert isinstance(result, GVAExtractionResult)
        assert (
            abs(result.g_V - true_gV) < 0.15
        ), f"g_V: expected {true_gV}, got {result.g_V}"

    def test_gV_gA_pure_gt(self):
        """For pure Gamow-Teller, g_V ≈ 0."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])

        # Pure GT: C(W) = g_A² (constant)
        true_gA = 0.8
        cw = true_gA**2 + np.random.normal(0, 0.02, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.02
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.02,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.02,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_gV_gA(cw_result, M_F=0.0, M_GT=1.0)

        assert isinstance(result, GVAExtractionResult)
        assert (
            abs(result.g_A - true_gA) < 0.15
        ), f"g_A: expected {true_gA}, got {result.g_A}"

    def test_gV_gA_mixed_transition(self):
        """For mixed transition, extract both g_V and g_A."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])

        # Mixed: C(W) = g_V² + g_A² (constant for simplicity)
        true_gV, true_gA = 0.4, 0.6
        cw = true_gV**2 + true_gA**2 + np.random.normal(0, 0.02, len(W))

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )

        extractor = CWExtractor(
            config, cw, energies_keV, flux_errors=np.ones_like(cw) * 0.02
        )
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.02,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.02,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_gV_gA(cw_result, M_F=1.0, M_GT=1.0)

        assert isinstance(result, GVAExtractionResult)
        # Check that sum of squares is correct
        reconstructed = result.g_V**2 + result.g_A**2
        expected = true_gV**2 + true_gA**2
        assert (
            abs(reconstructed - expected) < 0.15
        ), f"g_V²+g_A²: expected {expected}, got {reconstructed}"

    def test_gV_gA_summary(self):
        """summary() should return a non-empty string with key info."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw = np.full_like(W, 0.5)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )

        extractor = CWExtractor(config, cw, energies_keV)
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.02,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.02,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        result = extractor.fit_gV_gA(cw_result, M_F=1.0, M_GT=0.0)

        summary = result.summary()
        assert "g_V" in summary
        assert "g_A" in summary
        assert "χ²" in summary
        assert "p-value" in summary


class TestKurieAnalysis:
    """Test Kurie plot analysis."""

    def test_kurie_analysis_allowed_decay(self):
        """Kurie analysis on allowed decay should produce reasonable C(W)."""
        np.random.seed(42)
        energies_keV = np.linspace(5, 280, 80)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F1",
            use_screening=True,
            use_exchange=True,
            use_finite_size=True,
            use_radiative=True,
        )
        spectrum = BetaSpectrum.from_config(config)
        W, E_MeV = spectrum.get_energy_grid(config)
        flux_theory = spectrum(W)

        # Interpolate to our energy grid
        from scipy.interpolate import interp1d

        flux_interp = interp1d(
            E_MeV * 1000, flux_theory, kind="linear", fill_value=0.0, bounds_error=False
        )
        flux = flux_interp(energies_keV)

        # Add noise
        flux += np.random.normal(0, 0.02 * np.max(flux), len(flux))

        extractor = CWExtractor(config, flux, energies_keV)
        cw_result, fit_result = extractor.kurie_analysis(endpoint_keV=294.0)

        assert isinstance(cw_result, CWExtractionResult)
        assert cw_result.fit_result is not None
        assert np.all(cw_result.cw_values > 0)

    def test_kurie_plot_values(self):
        """Kurie plot should produce valid values."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 50)
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw = np.ones_like(W) * 1.0

        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        E_plot, kurie = cw_result.kurie_plot()

        assert len(E_plot) == len(energies_keV)
        assert len(kurie) == len(energies_keV)
        assert np.all(kurie >= 0), "Kurie values should be non-negative"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_not_enough_points_for_fit(self):
        """Fitting with too few points should raise ValueError."""
        np.random.seed(42)
        energies_keV = np.array([10.0, 20.0, 30.0])
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw = np.ones_like(W)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )

        extractor = CWExtractor(config, cw, energies_keV)
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        with pytest.raises(ValueError, match="Not enough data points"):
            extractor.fit_CW_parametrization(cw_result, parametrization="quadratic")

    def test_not_enough_points_for_gV_gA(self):
        """g_V/g_A fit with too few points should raise ValueError."""
        np.random.seed(42)
        energies_keV = np.array([10.0, 20.0, 30.0])
        W = np.array([T_to_W(e / 1000.0) for e in energies_keV])
        cw = np.ones_like(W)

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )

        extractor = CWExtractor(config, cw, energies_keV)
        cw_result = CWExtractionResult(
            energies_W=W.copy(),
            energies_keV=energies_keV.copy(),
            cw_values=cw.copy(),
            cw_errors=np.ones_like(cw) * 0.05,
            flux=cw.copy(),
            flux_errors=np.ones_like(cw) * 0.05,
            endpoint_W=T_to_W(0.294),
            endpoint_keV=294.0,
        )

        with pytest.raises(ValueError, match="Not enough data points"):
            extractor.fit_gV_gA(cw_result, M_F=1.0, M_GT=1.0)

    def test_default_uncertainties(self):
        """Default uncertainties should be sqrt(flux)."""
        np.random.seed(42)
        energies_keV = np.linspace(10, 250, 30)
        flux = np.array(
            [
                100.0,
                200.0,
                150.0,
                180.0,
                120.0,
                90.0,
                80.0,
                70.0,
                60.0,
                50.0,
                45.0,
                40.0,
                35.0,
                30.0,
                25.0,
                20.0,
                18.0,
                15.0,
                12.0,
                10.0,
                8.0,
                7.0,
                6.0,
                5.0,
                4.0,
                3.0,
                2.0,
                1.5,
                1.0,
                0.5,
            ]
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )

        extractor = CWExtractor(config, flux, energies_keV)
        expected_unc = np.sqrt(np.maximum(flux, 1.0))
        np.testing.assert_allclose(extractor.flux_errors, expected_unc)
