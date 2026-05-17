"""Tests for fitter.extractor."""

import numpy as np

from fitter.extractor import (
    CWExtractor,
    GVAExtractor,
    CWExtractionResult,
    GVAExtractionResult,
)


class TestCWExtractor:
    """Test CWExtractor."""

    def test_extract(self):
        from beta_spectrum import SpectrumConfig

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )
        extractor = CWExtractor(config)

        energies_keV = np.linspace(0.01, 0.29, 50)
        measured = np.linspace(10.0, 0.1, 50)
        errors = np.sqrt(measured)
        model = np.linspace(9.5, 0.09, 50)

        result = extractor.extract(
            measured_counts=measured,
            measured_errors=errors,
            energies_keV=energies_keV,
            model_counts=model,
            endpoint_keV=294.0,
        )

        assert isinstance(result, CWExtractionResult)
        assert len(result.cw_values) == 50
        assert len(result.energies_keV) == 50
        assert np.isclose(result.endpoint_keV, 294.0)
        assert np.all(np.isfinite(result.cw_values))

    def test_extract_with_zero_model(self):
        from beta_spectrum import SpectrumConfig

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )
        extractor = CWExtractor(config)

        energies_keV = np.linspace(0.01, 0.29, 50)
        measured = np.ones(50) * 5.0
        errors = np.ones(50) * 2.0
        model = np.zeros(50)  # All zeros

        result = extractor.extract(
            measured_counts=measured,
            measured_errors=errors,
            energies_keV=energies_keV,
            model_counts=model,
            endpoint_keV=294.0,
        )

        # C(W) should be 1.0 where model is zero (default)
        assert np.all(result.cw_values == 1.0)


class TestGVAExtractor:
    """Test GVAExtractor."""

    def test_extract_from_curve_basic(self):
        energies_W = np.linspace(1.1, 1.57, 100)
        endpoint_W = 1.575  # ~294 keV
        cw_values = np.ones_like(energies_W)  # C(W) = 1 for allowed

        result = GVAExtractor.extract_from_curve(
            cw_values=cw_values,
            energies_W=energies_W,
            endpoint_W=endpoint_W,
        )

        assert isinstance(result, GVAExtractionResult)
        assert np.isclose(result.gv_eff, 1.0)
        assert result.chi2_per_dof >= 0

    def test_extract_from_curve_with_slope(self):
        energies_W = np.linspace(1.1, 1.57, 100)
        endpoint_W = 1.575
        # C(W) with slight slope
        cw_values = 1.0 + 0.1 * (energies_W - endpoint_W)

        result = GVAExtractor.extract_from_curve(
            cw_values=cw_values,
            energies_W=energies_W,
            endpoint_W=endpoint_W,
        )

        assert isinstance(result, GVAExtractionResult)
        # g_A should be non-zero due to slope
        assert abs(result.ga_eff) > 0

    def test_extract_from_curve_insufficient_data(self):
        energies_W = np.array([1.0, 1.5])
        endpoint_W = 1.575
        cw_values = np.array([1.0, 1.0])

        result = GVAExtractor.extract_from_curve(
            cw_values=cw_values,
            energies_W=energies_W,
            endpoint_W=endpoint_W,
        )

        assert isinstance(result, GVAExtractionResult)
        assert result.chi2_per_dof == 0.0

    def test_gva_summary(self):
        result = GVAExtractionResult(
            gv_eff=1.0,
            gv_error=0.01,
            ga_eff=0.5,
            ga_error=0.02,
            correlation=-0.3,
            chi2_per_dof=1.2,
        )
        summary = result.summary()
        assert "g_V" in summary
        assert "g_A" in summary
        assert "1.0" in summary
        assert "0.5" in summary
        assert "1.2" in summary


class TestCWExtractionResult:
    """Test CWExtractionResult methods."""

    def test_kurie_plot(self):
        energies_W = np.linspace(1.1, 1.57, 100)
        endpoint_W = 1.575

        result = CWExtractionResult(
            energies_W=energies_W,
            energies_keV=energies_W * 510.998950 - 510.998950,
            cw_values=np.ones_like(energies_W),
            cw_errors=np.zeros_like(energies_W),
            flux=np.ones_like(energies_W) * 10.0,
            flux_errors=np.ones_like(energies_W) * 3.0,
            endpoint_W=endpoint_W,
            endpoint_keV=294.0,
        )

        E_keV, kurie = result.kurie_plot()
        assert len(E_keV) > 0
        assert len(kurie) > 0
        assert np.all(np.isfinite(kurie))
        assert np.all(kurie >= 0)
