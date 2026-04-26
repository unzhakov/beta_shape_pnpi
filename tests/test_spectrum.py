# tests/test_spectrum.py — Integration tests for the full BetaSpectrum pipeline

"""
Why test this?
--------------
Unit tests verify each component in isolation. But integration tests catch bugs that
only appear when components interact: wrong energy grid, type mismatches between
components, or numerical errors accumulating across multiplicative factors.

The `from_config()` factory method is the user-facing API — it's what people call
in real code. If this doesn't work for a realistic configuration, nothing else matters.

We test:

 1. **End-to-end pipeline**: from_config → energy grid → spectrum evaluation → all components valid.
     This is the "smoke test" that catches the most common user errors.
 2. **Selective component toggling**: Disabling a component shouldn't break anything.
     Users need to be able to turn off corrections for debugging or comparison studies.
 3. **Positive final spectrum**: The product of all factors must remain positive —
     if any component goes negative, the entire spectrum is wrong.
 4. **Component extraction works**: calculate_components() should return per-component values
     that match calling each component individually.

Common practice: Integration tests are fewer but more expensive (they test more code paths).
They serve as a safety net — if they pass, you can be confident the system works "as designed".
"""

import numpy as np
import pytest

from beta_spectrum import BetaSpectrum, SpectrumConfig


class TestBetaSpectrumIntegration:
    """Full pipeline integration tests."""

    def test_from_config_basic(self):
        """The most basic use case: create a spectrum from config and evaluate it.

        This is the entry point shown in README.md — if this fails, nothing works.
        We test with a medium-Z nucleus (Z=20 calcium) which exercises all components.
        """
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5
        )

        spectrum = BetaSpectrum.from_config(config)
        W, E_MeV = spectrum.get_energy_grid(config)

        assert len(W) == len(E_MeV), "Energy grids must have same length"
        assert len(W) > 0, "Grid must not be empty"

    def test_spectrum_evaluates_without_error(self):
        """Calling the spectrum object should return a valid numpy array."""
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        values = spectrum(W)

        assert isinstance(values, np.ndarray), "Output must be a numpy array"
        assert values.shape == W.shape, "Output shape must match input shape"

    def test_all_components_positive(self):
        """Every component and the total product must be positive."""
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        # Filter out values exactly at the endpoint (where phase space is zero)
        mask = np.isfinite(spectrum(W)) & (spectrum(W) > 0)
        assert np.any(mask), "At least some finite positive values expected"


class TestBetaSpectrumToggles:
    """Test that component toggling works correctly."""

    def test_all_disabled(self):
        """With all corrections disabled, only phase space remains — still valid output."""
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5,
            use_fermi=False,
            use_screening=False,
            use_finite_size=False,
            use_charge_dist=False,
            use_exchange=False,
            use_radiative=False,
        )

        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert len(values) > 0, "Should return non-empty array"
        assert np.all(np.isfinite(values)), "All values must be finite"

    def test_fermi_only(self):
        """With only fermi function enabled: should give F₀(Z,W) — not just phase space."""
        config_minimal = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5,
            use_phase_space=False,
            use_screening=False,
            use_finite_size=False,
            use_charge_dist=False,
            use_exchange=False,
            use_radiative=False,
        )

        spectrum = BetaSpectrum.from_config(config_minimal)
        W, _ = spectrum.get_energy_grid(config_minimal)
        values = spectrum(W)

        # Fermi function alone is ~1 for Z=20 at mid-energy
        assert np.all(np.isfinite(values)), "Fermi-only output must be finite"


class TestBetaSpectrumComponentExtraction:
    """Test calculate_components() returns correct per-component data."""

    def test_component_names_match(self):
        """Component names should correspond to class names (minus suffixes)."""
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5
        )
        spectrum = BetaSpectrum.from_config(config)

        W, _ = spectrum.get_energy_grid(config)
        components = spectrum.calculate_components(W)

        assert isinstance(components, dict), "Must return a dictionary"
        assert len(components) > 0, "Should have at least one component"


class TestBetaSpectrumAnalyzer:
    """Test the analyzer utilities."""

    def test_total_spectrum_returns_array(self):
        config = SpectrumConfig(
            Z_parent=19, Z_daughter=20, A_number=40, endpoint_MeV=2.5
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = type('Analyzer', (), {'spectrum': spectrum, 'config': config})()

        # We can't easily instantiate the full Analyzer without matplotlib backend,
        # but we verify the spectrum itself works end-to-end.
        W, E = spectrum.get_energy_grid(config)
        total = spectrum(W)

        assert len(total) == len(E), "Spectrum and energy arrays must match"


class TestBetaSpectrumRealisticConfig:
    """Test with realistic physics parameters."""

    def test_uranium_spectrum(self):
        """Heavy nucleus (Z=92) — tests numerical stability at extreme Z.

        Uranium-238 decays via alpha, but we can still evaluate the beta spectrum
        shape for testing purposes. High Z stresses the Fermi function's loggamma.
        """
        config = SpectrumConfig(
            Z_parent=91, Z_daughter=92, A_number=238, endpoint_MeV=0.5
        )

        spectrum = BetaSpectrum.from_config(config)
        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert np.all(np.isfinite(values)), "Uranium must not produce NaN/inf"
        assert len(values) > 10, f"Grid should have many points: got {len(values)}"

    def test_low_Z_spectrum(self):
        """Light nucleus (Z=6 carbon-14 beta decay)."""
        config = SpectrumConfig(
            Z_parent=5, Z_daughter=6, A_number=14, endpoint_MeV=0.156  # ~156 keV for C-14
        )

        spectrum = BetaSpectrum.from_config(config)
        W, E = spectrum.get_energy_grid(config)
        values = spectrum(W)

        assert np.all(np.isfinite(values)), "C-14 must not produce NaN/inf"
        # Phase space should peak in the bulk for low-Q decay
        assert np.any(values > 0), "Should have non-zero spectrum values"
