"""
Physics tests for RadiativeCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.radiative import RadiativeCorrection


class TestRadiativeEndpointHandling:
    """Test how endpoint divergence is handled."""

    @pytest.mark.physics
    def test_resummation_finite_at_endpoint(self):
        """With resummation enabled, R(W₀−ε) must be finite for small ε."""
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        W_near_end = np.array([4.999])
        result = rc(W_near_end)

        assert not np.isnan(result[0]), "Resummed correction must be finite at endpoint"
        assert not np.isinf(result[0]), "Resummed correction must not diverge at endpoint"

    @pytest.mark.physics
    def test_standard_mode_diverges(self):
        """Without resummation, δᵣ → −∞ as W → W₀ (logarithmic divergence)."""
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        W_far = np.array([4.9])   # ΔW = 0.1
        W_close = np.array([4.99])  # ΔW = 0.01

        val_far = rc(W_far)[0]
        val_close = rc(W_close)[0]

        assert val_close < val_far, "Standard mode must show divergent behavior near endpoint"


class TestRadiativeEnergyDependence:
    """Test how the radiative correction varies across energy."""

    @pytest.mark.physics
    def test_correction_grows_toward_endpoint(self):
        """δᵣ grows as W → W₀ due to ln(W₀−W) divergence."""
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)

        mid_W = np.array([3.0])
        near_end_W = np.array([4.8])

        r_mid = rc(mid_W)[0]
        r_near = rc(near_end_W)[0]

        assert r_near <= r_mid, f"Correction factor should decrease near endpoint: {r_mid} → {r_near}"


class TestRadiativeResummationSwitch:
    """Test that the resummation flag actually changes behavior."""

    @pytest.mark.physics
    def test_resummed_differs_from_standard(self):
        """Near endpoint (within delta_cut), resummed and standard modes should differ."""
        rc_resummed = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        rc_standard = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        # Far from endpoint: both identical
        W_far = np.array([3.0])
        assert np.isclose(rc_resummed(W_far)[0], rc_standard(W_far)[0], rtol=1e-10)

        # Near endpoint: resummed stays finite, standard diverges
        W_near = np.array([4.9995])
        r_resummed = rc_resummed(W_near)[0]
        r_standard = rc_standard(W_near)[0]

        assert not np.isnan(r_resummed), "Resummed mode must be finite"
        assert not np.isnan(r_standard), "Standard mode must be finite (mask protects endpoint)"
        assert r_resummed > r_standard, f"Near endpoint: resummed={r_resummed} > standard={r_standard}"

    @pytest.mark.physics
    def test_resummation_only_near_endpoint(self):
        """Per spec, resummation only applies when (W0 - W) < delta_cut."""
        rc_resummed = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        rc_standard = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        # delta_cut defaults to 1e-3, so W0 - W = 0.01 > delta_cut
        W_far = np.array([4.99])
        r_resummed_far = rc_resummed(W_far)[0]
        r_standard_far = rc_standard(W_far)[0]

        assert np.isclose(r_resummed_far, r_standard_far, rtol=1e-10), (
            f"Far from endpoint: resummed={r_resummed_far}, standard={r_standard_far}"
        )

        # Very near endpoint (delta_W < delta_cut): resummed is finite and larger
        W_near = np.array([4.9995])
        r_resummed_near = rc_resummed(W_near)[0]
        r_standard_near = rc_standard(W_near)[0]

        assert not np.isnan(r_resummed_near)
        assert not np.isnan(r_standard_near)
        assert r_resummed_near > r_standard_near


class TestRadiativeZDependence:
    """Test Z-dependent O(Z*alpha^2) correction from Sirlin 1987."""

    @pytest.mark.physics
    def test_correction_increases_with_Z(self):
        """Higher Z produces larger O(Z*alpha^2) correction."""
        rc_low = RadiativeCorrection(W0=5.0, Z=20, use_endpoint_resummation=True)
        rc_high = RadiativeCorrection(W0=5.0, Z=80, use_endpoint_resummation=True)

        W_test = np.array([3.0])
        assert rc_high(W_test)[0] > rc_low(W_test)[0]

    @pytest.mark.physics
    def test_different_Z_different_results(self):
        """Different Z values must produce measurably different corrections."""
        rc_z1 = RadiativeCorrection(W0=5.0, Z=1, use_endpoint_resummation=True)
        rc_z50 = RadiativeCorrection(W0=5.0, Z=50, use_endpoint_resummation=True)

        W_test = np.linspace(1.5, 4.5, 10)
        diff = np.abs(rc_z50(W_test) - rc_z1(W_test))
        assert np.max(diff) > 1e-4, f"Z=1 and Z=50 should differ: max diff = {np.max(diff)}"

    @pytest.mark.physics
    def test_z_correction_magnitude_reasonable(self):
        """For Z=92 (uranium), O(Z*alpha^2) should be reasonable."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)
        W_test = np.array([3.0])
        result = rc(W_test)
        assert 0.9 < result[0] < 3.0, f"Z=92 correction should be reasonable: {result[0]}"


class TestRadiativeAParameter:
    """Test A (mass number) parameter for nuclear model."""

    @pytest.mark.physics
    def test_a_affects_nuclear_model_correction(self):
        """A parameter affects the nuclear-structure-dependent part of O(Z*alpha^2)."""
        rc_a100 = RadiativeCorrection(W0=5.0, Z=50, A=100, use_endpoint_resummation=True)
        rc_a140 = RadiativeCorrection(W0=5.0, Z=50, A=140, use_endpoint_resummation=True)

        W_test = np.array([3.0])
        assert rc_a100(W_test)[0] != rc_a140(W_test)[0], (
            f"Different A values must differ: A=100 -> {rc_a100(W_test)[0]}, "
            f"A=140 -> {rc_a140(W_test)[0]}"
        )
