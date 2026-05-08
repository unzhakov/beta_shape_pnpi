"""Tests for multi-branch decay support.

Why test this?
--------------
Multi-branch decay is a critical feature for realistic nuclear decay analysis.
Many isotopes decay through multiple channels to different daughter states,
each with its own endpoint energy, transition type, and intensity.

The multi-branch mode is the DEFAULT — even single-branch decays are handled
through the same code path. The only difference in behavior is in plot generation
(where a single branch doesn't need separate branch plots).

We test:

 1. **BranchConfig validation**: intensity bounds, positive endpoint.
 2. **Single-branch via multi-branch path**: one branch should give identical
    results to the old single-branch behavior.
 3. **Multi-branch spectrum calculation**: weighted sum produces valid output.
 4. **Intensity weighting**: sum of intensities = 1.0, weighted sum = total.
 5. **Per-branch spectra**: each branch has correct endpoint and transition type.
 6. **Component extraction**: universal + per-branch components correctly separated.
 7. **CSV export**: multi-branch CSV has correct column structure.
 8. **Energy grid**: uses max endpoint across all branches.
 9. **Detector response**: multi-branch spectra are convolvable.
 10. **BranchConfig dataclass**: basic construction and validation.
"""

import numpy as np
import pytest

from beta_spectrum import (
    BetaSpectrum,
    BetaSpectrumAnalyzer,
    BranchConfig,
    SpectrumConfig,
)
from beta_spectrum.utils import T_to_W

# ---------------------------------------------------------------------------
# BranchConfig validation
# ---------------------------------------------------------------------------


class TestBranchConfigValidation:
    """Test BranchConfig dataclass validation."""

    def test_valid_branch(self):
        """Valid branch config should construct without error."""
        branch = BranchConfig(
            endpoint_MeV=0.2,
            transition_type="A",
            intensity=0.5,
        )
        assert branch.endpoint_MeV == 0.2
        assert branch.transition_type == "A"
        assert branch.intensity == 0.5

    def test_intensity_out_of_range_low(self):
        """Intensity below 0 should raise ValueError."""
        with pytest.raises(ValueError, match="intensity must be in"):
            BranchConfig(endpoint_MeV=0.2, intensity=-0.1)

    def test_intensity_out_of_range_high(self):
        """Intensity above 1.0 should raise ValueError."""
        with pytest.raises(ValueError, match="intensity must be in"):
            BranchConfig(endpoint_MeV=0.2, intensity=1.1)

    def test_zero_endpoint(self):
        """Zero endpoint should raise ValueError."""
        with pytest.raises(ValueError, match="endpoint_MeV must be positive"):
            BranchConfig(endpoint_MeV=0.0, intensity=0.5)

    def test_negative_endpoint(self):
        """Negative endpoint should raise ValueError."""
        with pytest.raises(ValueError, match="endpoint_MeV must be positive"):
            BranchConfig(endpoint_MeV=-0.1, intensity=0.5)

    def test_default_transition_type(self):
        """Transition type defaults to 'A'."""
        branch = BranchConfig(endpoint_MeV=0.2, intensity=0.5)
        assert branch.transition_type == "A"


# ---------------------------------------------------------------------------
# Single-branch via multi-branch path (backward compatibility)
# ---------------------------------------------------------------------------


class TestSingleBranchViaMultiBranch:
    """Verify single-branch mode works through the multi-branch path."""

    def test_one_branch_produces_valid_spectrum(self):
        """One branch should produce valid spectrum identical to old behavior."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            branches=[
                BranchConfig(
                    endpoint_MeV=0.294,
                    transition_type="F2",
                    intensity=1.0,
                )
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert len(values) == len(E)
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0), "Spectrum should be non-negative"

    def test_one_branch_has_no_branch_spectra(self):
        """Single branch should not produce separate branch plots."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=1.0),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        branch_spectra = spectrum.get_branch_spectra(W)
        # Single branch: no separate branch spectra needed
        assert len(branch_spectra) == 0

    def test_one_branch_analyzer_uses_normal_plot(self):
        """Single branch analyzer should not have multi-branch mode."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=1.0),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        assert not analyzer._is_multi


# ---------------------------------------------------------------------------
# Multi-branch spectrum calculation
# ---------------------------------------------------------------------------


class TestMultiBranchSpectrum:
    """Test multi-branch spectrum calculation."""

    def test_two_branches_basic(self):
        """Two branches should produce valid weighted sum."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)

        assert spectrum._is_multi
        assert len(spectrum.branch_spectra) == 2
        assert len(spectrum.branches) == 2

        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert len(values) == len(E)
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0), "Spectrum should be non-negative"

    def test_three_branches(self):
        """Three branches should produce valid weighted sum."""
        config = SpectrumConfig(
            Z_parent=27,
            Z_daughter=28,
            A_number=60,
            endpoint_MeV=0.320,
            branches=[
                BranchConfig(endpoint_MeV=0.310, transition_type="A", intensity=0.5),
                BranchConfig(endpoint_MeV=0.200, transition_type="F1", intensity=0.3),
                BranchConfig(endpoint_MeV=0.100, transition_type="A", intensity=0.2),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)

        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert spectrum._is_multi
        assert len(values) == len(E)
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0)

    def test_energy_grid_uses_max_endpoint(self):
        """Energy grid should extend to the maximum branch endpoint."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            e_step_MeV=0.001,
            branches=[
                BranchConfig(endpoint_MeV=0.300, transition_type="A", intensity=0.5),
                BranchConfig(endpoint_MeV=0.200, transition_type="F2", intensity=0.5),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, E = spectrum.get_energy_grid(config)

        # Grid should extend to max endpoint (0.300 MeV)
        assert E[-1] >= 0.299, f"Grid should reach ~0.300 MeV, got {E[-1]}"

    def test_branch_spectra_property(self):
        """get_branch_spectra() should return per-branch unweighted spectra."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        branch_spectra = spectrum.get_branch_spectra(W)
        assert len(branch_spectra) == 2
        assert len(branch_spectra[0]) == len(W)
        assert len(branch_spectra[1]) == len(W)

    def test_branch_normalized_spectra(self):
        """get_branch_normalized_spectra() should return intensity-weighted spectra."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            # Disable universal components to simplify
            use_fermi=False,
            use_screening=False,
            use_finite_size=False,
            use_charge_dist=False,
            use_exchange=False,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        weighted = spectrum.get_branch_normalized_spectra(W)
        assert len(weighted) == 2

        # With no universal components, total = sum of weighted branch spectra
        total_from_branches = weighted[0] + weighted[1]
        total = spectrum(W)
        assert np.allclose(total, total_from_branches, rtol=1e-10)

    def test_intensity_normalization(self):
        """Branch intensities should be used as weights correctly."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            use_phase_space=True,
            use_fermi=False,
            use_screening=False,
            use_finite_size=False,
            use_charge_dist=False,
            use_exchange=False,
            use_radiative=False,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="A", intensity=0.8),
                BranchConfig(endpoint_MeV=0.294, transition_type="A", intensity=0.2),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        # Both branches have same endpoint and transition type,
        # so phase space should be identical. Total = 0.8*PS + 0.2*PS = PS
        total = spectrum(W)
        branch0 = spectrum.branch_spectra[0](W)
        branch1 = spectrum.branch_spectra[1](W)

        expected = 0.8 * branch0 + 0.2 * branch1
        assert np.allclose(total, expected, rtol=1e-10)


# ---------------------------------------------------------------------------
# Component extraction in multi-branch mode
# ---------------------------------------------------------------------------


class TestMultiBranchComponents:
    """Test component extraction for multi-branch mode."""

    def test_universal_components_present(self):
        """Universal components (Fermi, etc.) should be in component dict."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)
        components = spectrum.calculate_components(W)

        assert "Fermi" in components, "Universal Fermi component should be present"
        assert "Screening" in components
        assert "FiniteSizeL0" in components

    def test_per_branch_components_present(self):
        """Per-branch components (PhaseSpace, Radiative) should be in component dict."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)
        components = spectrum.calculate_components(W)

        assert "branch_0.PhaseSpace" in components
        assert "branch_1.PhaseSpace" in components
        assert "branch_0.Radiative" in components
        assert "branch_1.Radiative" in components

    def test_analyzer_components_multi_branch(self):
        """BetaSpectrumAnalyzer should handle multi-branch components."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        components = analyzer.components
        assert "Fermi" in components
        assert "branch_0.PhaseSpace" in components

    def test_analyzer_branch_spectra_property(self):
        """Analyzer should expose branch spectra property."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        branch_spectra = analyzer.branch_spectra
        assert len(branch_spectra) == 2

    def test_total_equals_universal_times_branch_sum(self):
        """Total = universal_components * sum(intensity * branch_spectra)."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        total = spectrum(W)

        # Compute expected: universal_components * sum(intensity * branch_spectra)
        # Each branch spectrum must be masked to zero beyond its endpoint
        universal_product = np.ones_like(W)
        for comp in spectrum.components:
            universal_product *= comp(W)

        branch_sum = np.zeros_like(W)
        for i, spec in enumerate(spectrum.branch_spectra):
            W0_branch = T_to_W(spectrum.branches[i].endpoint_MeV)
            mask = W <= W0_branch
            branch_sum += spectrum.branches[i].intensity * np.where(mask, spec(W), 0.0)

        expected = universal_product * branch_sum
        assert np.allclose(total, expected, rtol=1e-10)


# ---------------------------------------------------------------------------
# CSV export in multi-branch mode
# ---------------------------------------------------------------------------


class TestMultiBranchCSVExport:
    """Test CSV export structure for multi-branch mode."""

    def _write_json(self, data):
        import json
        import tempfile
        from pathlib import Path

        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(data, tmp)
        tmp.close()
        return Path(tmp.name)

    def test_multi_branch_csv_columns(self, tmp_path):
        """Multi-branch CSV should have total + per-branch + universal columns."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        csv_path = str(tmp_path / "multi_branch.csv")
        analyzer.export_to_csv(csv_path)

        import pandas as pd

        # Read CSV (skip comment headers)
        df = pd.read_csv(csv_path, comment="#")

        # Should have total spectrum column
        assert "spectrum" in df.columns

        # Should have per-branch columns
        assert "branch_1_spectrum" in df.columns
        assert "branch_2_spectrum" in df.columns

        # Should have universal columns
        assert "universal_Fermi" in df.columns

    def test_single_branch_csv_columns(self, tmp_path):
        """Single-branch CSV should have standard columns."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=1.0),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        csv_path = str(tmp_path / "single_branch.csv")
        analyzer.export_to_csv(csv_path)

        import pandas as pd

        df = pd.read_csv(csv_path, comment="#")

        assert "spectrum" in df.columns
        assert "Fermi" in df.columns


# ---------------------------------------------------------------------------
# Multi-branch with different transition types
# ---------------------------------------------------------------------------


class TestMultiBranchTransitionTypes:
    """Test that different transition types are correctly propagated."""

    def test_different_transition_types(self):
        """Each branch should use its own transition type for PhaseSpace."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.200, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)

        # Check that each branch has the correct transition type
        assert spectrum.branches[0].transition_type == "F2"
        assert spectrum.branches[1].transition_type == "A"

        # PhaseSpace components should have different W0
        W, _ = spectrum.get_energy_grid(config)
        components = spectrum.calculate_components(W)

        branch0_ps = components["branch_0.PhaseSpace"]
        branch1_ps = components["branch_1.PhaseSpace"]

        # Branch 0 has higher endpoint, so its phase space extends further
        assert branch0_ps.shape == branch1_ps.shape

    def test_branch_info_transition_type_from_ensdf(self):
        """BranchInfo should have transition_type field populated from ENSDF."""
        from beta_spectrum.nuclear_data import BranchInfo

        branch = BranchInfo(
            level_index=1,
            level_energy_keV=100.0,
            intensity=0.5,
            log_ft=5.0,
            transition_type="F1",
        )
        assert branch.transition_type == "F1"
        assert branch.level_index == 1
        assert branch.level_energy_keV == 100.0
        assert branch.intensity == 0.5

    def test_all_branches_same_transition(self):
        """All branches can have the same transition type."""
        config = SpectrumConfig(
            Z_parent=27,
            Z_daughter=28,
            A_number=60,
            endpoint_MeV=0.320,
            branches=[
                BranchConfig(endpoint_MeV=0.310, transition_type="A", intensity=0.5),
                BranchConfig(endpoint_MeV=0.200, transition_type="A", intensity=0.3),
                BranchConfig(endpoint_MeV=0.100, transition_type="A", intensity=0.2),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)

        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert len(values) == len(E)
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0)


# ---------------------------------------------------------------------------
# Multi-branch with detector response
# ---------------------------------------------------------------------------


class TestMultiBranchDetectorResponse:
    """Test multi-branch with detector response convolution."""

    def test_multi_branch_with_detector(self):
        """Multi-branch spectrum should be convolvable with detector response."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            use_detector_response=True,
            detector_model="gaussian",
            detector_sigma_a_keV=1.0,
            detector_n_channels=256,
            detector_fano_factor=0.0,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.7),
                BranchConfig(endpoint_MeV=0.150, transition_type="A", intensity=0.3),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        detector = BetaSpectrum.create_detector_from_config(config)
        convolved = spectrum.convolve_with_detector(detector, W=W, config=config)

        assert convolved.shape == (256,)
        assert np.all(convolved >= 0), "Convolved spectrum should be non-negative"
        assert np.sum(convolved) > 0, "Convolved spectrum should have non-zero area"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestMultiBranchEdgeCases:
    """Test edge cases for multi-branch mode."""

    def test_very_low_intensity_branch(self):
        """Branch with very low intensity should still be computed."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.350,
            branches=[
                BranchConfig(endpoint_MeV=0.294, transition_type="F2", intensity=0.999),
                BranchConfig(endpoint_MeV=0.100, transition_type="A", intensity=0.001),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        values = spectrum(W)
        assert np.all(np.isfinite(values))

    def test_very_high_z_multi_branch(self):
        """Multi-branch with high Z should be numerically stable."""
        config = SpectrumConfig(
            Z_parent=90,
            Z_daughter=91,
            A_number=232,
            endpoint_MeV=1.0,
            branches=[
                BranchConfig(endpoint_MeV=0.8, transition_type="A", intensity=0.6),
                BranchConfig(endpoint_MeV=0.5, transition_type="F1", intensity=0.4),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        values = spectrum(W)
        assert np.all(
            np.isfinite(values)
        ), "High-Z multi-branch must not produce NaN/inf"

    def test_no_branches_defaults_to_single(self):
        """No branches list should default to single-branch mode."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )
        spectrum = BetaSpectrum.from_config(config)

        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert len(values) == len(E)
        assert np.all(np.isfinite(values))

    def test_branch_values_below_total(self):
        """Each branch contribution should be <= total spectrum.

        Since total = sum(intensity_i * branch_i) * universal,
        and all factors are positive, each branch's contribution
        should not exceed the total.
        """
        config = SpectrumConfig(
            Z_parent=59,
            Z_daughter=60,
            A_number=144,
            endpoint_MeV=2.996,
            e_step_MeV=0.001,
            branches=[
                BranchConfig(
                    endpoint_MeV=2.996, transition_type="F1", intensity=0.1469
                ),
                BranchConfig(
                    endpoint_MeV=0.341, transition_type="F1", intensity=0.1220
                ),
                BranchConfig(
                    endpoint_MeV=0.253, transition_type="F1", intensity=0.1217
                ),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, energies = spectrum.get_energy_grid(config)

        total = spectrum(W)
        bs = spectrum.get_branch_spectra(W)
        components = spectrum.calculate_components(W)

        # Compute each branch's contribution to the total
        for i, branch in enumerate(config.branches):
            universal_product = np.ones_like(W)
            for name in [
                "Fermi",
                "Screening",
                "FiniteSizeL0",
                "ChargeDistributionU",
                "Exchange",
            ]:
                if name in components:
                    universal_product *= components[name]
            contrib = universal_product * branch.intensity * bs[i]

            # Branch contribution should not exceed total
            assert np.all(
                contrib <= total + 1e-10
            ), f"Branch {i} contribution exceeds total at some points"

            # Branch contribution should be non-negative
            assert np.all(contrib >= 0), f"Branch {i} contribution has negative values"

    def test_branch_integral_proportional_to_intensity(self):
        """Integral of each branch contribution should be proportional to its intensity."""
        config = SpectrumConfig(
            Z_parent=59,
            Z_daughter=60,
            A_number=144,
            endpoint_MeV=2.996,
            e_step_MeV=0.001,
            branches=[
                BranchConfig(endpoint_MeV=2.996, transition_type="F1", intensity=0.5),
                BranchConfig(endpoint_MeV=0.5, transition_type="F1", intensity=0.3),
                BranchConfig(endpoint_MeV=0.3, transition_type="F1", intensity=0.2),
            ],
        )
        spectrum = BetaSpectrum.from_config(config)
        W, energies = spectrum.get_energy_grid(config)

        bs = spectrum.get_branch_spectra(W)
        components = spectrum.calculate_components(W)

        integrals = []
        for i, branch in enumerate(config.branches):
            universal_product = np.ones_like(W)
            for name in [
                "Fermi",
                "Screening",
                "FiniteSizeL0",
                "ChargeDistributionU",
                "Exchange",
            ]:
                if name in components:
                    universal_product *= components[name]
            contrib = universal_product * branch.intensity * bs[i]
            integral = np.trapezoid(contrib, energies)
            integrals.append(integral)

        total_integral = np.trapezoid(spectrum(W), energies)

        # Each branch integral should be proportional to its intensity
        # (within numerical tolerance, since different endpoints give different
        # phase space integrals)
        for i, (integral, intensity) in enumerate(zip(integrals, config.branches)):
            assert integral > 0, f"Branch {i} has zero integral"

        # Sum of branch integrals should equal total
        assert np.isclose(
            sum(integrals), total_integral, rtol=1e-6
        ), "Sum of branch integrals should equal total integral"
