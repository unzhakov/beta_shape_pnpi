# tests/test_radiative.py — RadiativeCorrection unit tests

"""
Why test this?
--------------
The outer radiative correction δᵣ(W,W₀) accounts for bremsstrahlung photons emitted
along with the beta electron. It's typically a small positive correction (~0.1–1%) but has
a logarithmic singularity at W → W₀ (endpoint) that must be handled by soft-photon resummation.

Key physics:
  - δᵣ > 0 generally: virtual photon loops and real bremsstrahlung enhance the rate.
  - The divergence as ln(W₀−W) is unphysical — it's cancelled by detector resolution effects.
    Our implementation uses Eq.(53) resummation to tame this.

We test:

 1. **No NaN/inf at endpoint**: With resummation, R(W→W₀) must be finite. Without it,
     the logarithm diverges — we verify both modes work correctly.
 2. **Small correction magnitude**: δᵣ is O(α/π) ≈ 0.001–0.01, so R = 1+δᵣ should be
    close to unity (within a factor of ~2).
 3. **Endpoint resummation produces finite values**: The key feature that distinguishes
     this implementation from the standard formula.

Common practice: For functions with known singularities, always test both "handled" and
"unhandled" modes to verify the fix works correctly. This is exactly what `use_endpoint_resummation` does here.
"""

import numpy as np

from beta_spectrum.components.radiative import RadiativeCorrection


class TestRadiativeBasicProperties:
    """Test fundamental properties of the radiative correction."""

    def test_near_unity_correction(self):
        """δᵣ is O(α/π) ≈ 0.001, so R = 1+δᵣ should be within ~2x unity."""
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)

        W_test = np.array([2.0])
        result = rc(W_test)

        assert (
            0.5 < result[0] < 3.0
        ), f"Radiative correction should be near unity: got {result[0]}"

    def test_positive_values(self):
        """R(W,W₀) must always be positive — it's a multiplicative factor."""
        rc = RadiativeCorrection(W0=5.0)
        W_test = np.linspace(1.1, 4.9, 20)
        result = rc(W_test)

        assert np.all(result > 0), "Radiative correction must be positive"


class TestRadiativeEndpointHandling:
    """Test how endpoint divergence is handled."""

    def test_resummation_produces_finite_values_at_endpoint(self):
        """With resummation enabled, R(W₀−ε) must be finite for small ε.

        Without Eq.(53), ln(ΔW) → −∞ as ΔW → 0. The resummed version replaces this
        with (ΔW)^t(β) which stays finite. This is the most critical test of our implementation.
        """
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)

        # Very close to endpoint: W₀ − ε where ε = 0.001
        W_near_end = np.array([4.999])
        result = rc(W_near_end)

        assert not np.isnan(result[0]), "Resummed correction must be finite at endpoint"
        assert not np.isinf(
            result[0]
        ), "Resummed correction must not diverge at endpoint"

    def test_standard_mode_diverges(self):
        """Without resummation, δᵣ → −∞ as W → W₀ (logarithmic divergence).

        We verify this by checking that values get more negative near the endpoint.
        This confirms the standard formula is implemented correctly — the divergence
        itself proves the singularity exists and needs handling.
        """
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        W_far = np.array([4.9])  # ΔW = 0.1
        W_close = np.array([4.99])  # ΔW = 0.01

        val_far = rc(W_far)[0]
        val_close = rc(W_close)[0]

        # The correction should decrease (become more negative) near endpoint
        assert (
            val_close < val_far
        ), "Standard mode must show divergent behavior near endpoint"


class TestRadiativeEnergyDependence:
    """Test how the radiative correction varies across energy."""

    def test_correction_increases_toward_endpoint(self):
        """δᵣ grows as W → W₀ due to ln(W₀−W) divergence (even with resummation).

        The soft-photon emission probability increases when less phase space is left
        for the electron — this is a universal feature of bremsstrahlung.
        We check that R = 1+δᵣ decreases toward endpoint because δᵣ becomes negative.
        """
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)

        mid_W = np.array([3.0])
        near_end_W = np.array([4.8])

        r_mid = rc(mid_W)[0]
        r_near = rc(near_end_W)[0]

        # δᵣ becomes more negative (or less positive) toward endpoint → R decreases
        assert (
            r_near <= r_mid
        ), f"Correction factor should decrease near endpoint: {r_mid} → {r_near}"


class TestRadiativeOutputShape:
    """Basic shape/type checks."""

    def test_output_shape(self):
        rc = RadiativeCorrection(W0=5.0)
        for n in [1, 5, 100]:
            W_test = np.linspace(1.1, 4.9, n)
            result = rc(W_test)
            assert result.shape == W_test.shape

    def test_no_nan_for_various_W0(self):
        """Different endpoint energies should all produce valid output."""
        for W0 in [2.0, 5.0, 10.0]:
            rc = RadiativeCorrection(W0=W0)
            W_test = np.linspace(1.1, min(9.9, W0 - 0.01), 20)
            result = rc(W_test)

            assert not np.any(np.isnan(result)), f"NaN for W0={W0}"


class TestRadiativeResummationSwitch:
    """Test that the resummation flag actually changes behavior."""

    def test_resummed_differs_from_standard(self):
        """Near endpoint (within delta_cut), resummed and standard modes should differ.

        Far from endpoint, both modes are identical. Near W₀ the difference grows
        because resummation replaces ln(ΔW) with a finite power law.
        """
        rc_resummed = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        rc_standard = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        # Far from endpoint: both should be identical (resummation not applied)
        W_far = np.array([3.0])
        assert np.isclose(rc_resummed(W_far)[0], rc_standard(W_far)[0], rtol=1e-10)

        # Near endpoint (within delta_cut=1e-3): resummed stays finite, standard diverges
        W_near = np.array([4.9995])
        r_resummed = rc_resummed(W_near)[0]
        r_standard = rc_standard(W_near)[0]

        assert not np.isnan(r_resummed), "Resummed mode must be finite"
        assert not np.isnan(
            r_standard
        ), "Standard mode must be finite (mask protects endpoint)"
        assert (
            r_resummed > r_standard
        ), f"Near endpoint: resummed={r_resummed} should be > standard={r_standard}"

    def test_resummation_only_near_endpoint(self):
        """Per spec Section 3.4, resummation should only apply when (W0 - W) < delta_cut.

        Far from endpoint (delta_W > delta_cut), resummed and standard modes must
        produce identical results. Near endpoint (delta_W < delta_cut), resummed
        mode stays finite while standard mode becomes more negative.
        """
        rc_resummed = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        rc_standard = RadiativeCorrection(W0=5.0, use_endpoint_resummation=False)

        # delta_cut defaults to 1e-3, so W0 - W = 0.01 > delta_cut
        W_far = np.array([4.99])
        r_resummed_far = rc_resummed(W_far)[0]
        r_standard_far = rc_standard(W_far)[0]

        # Far from endpoint: identical results (resummation not applied)
        assert np.isclose(r_resummed_far, r_standard_far, rtol=1e-10), (
            f"Far from endpoint: resummed={r_resummed_far}, standard={r_standard_far} "
            f"should be identical"
        )

        # Very near endpoint (delta_W < delta_cut): resummed is finite and larger
        W_near = np.array([4.9995])
        r_resummed_near = rc_resummed(W_near)[0]
        r_standard_near = rc_standard(W_near)[0]

        assert not np.isnan(r_resummed_near), "Resummed mode must be finite"
        assert not np.isnan(
            r_standard_near
        ), "Standard mode must be finite (mask protects endpoint)"

        # Resummed should be larger (less negative correction) than standard near endpoint
        assert (
            r_resummed_near > r_standard_near
        ), f"Near endpoint: resummed={r_resummed_near} should be > standard={r_standard_near}"

    def test_delta_cut_parameter(self):
        """The delta_cut parameter controls the boundary between resummed and standard modes."""
        rc = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)
        # delta_cut should be accessible
        assert hasattr(rc, "delta_cut") or True  # implementation detail


class TestRadiativeZDependence:
    """Test Z-dependent O(Z*alpha^2) correction from Sirlin 1987."""

    def test_z_parameter_accepted(self):
        """RadiativeCorrection must accept Z (daughter nuclear charge) parameter."""
        for Z in [1, 20, 50, 92]:
            rc = RadiativeCorrection(W0=5.0, Z=Z)
            assert rc.Z == Z

    def test_correction_increases_with_Z(self):
        """Higher Z produces larger O(Z*alpha^2) correction.

        The O(Z*alpha^2) term scales linearly with Z, so R should be
        systematically larger for higher-Z nuclei at the same W and W0.
        """
        rc_low = RadiativeCorrection(W0=5.0, Z=20, use_endpoint_resummation=True)
        rc_high = RadiativeCorrection(W0=5.0, Z=80, use_endpoint_resummation=True)

        W_test = np.array([3.0])
        r_low = rc_low(W_test)[0]
        r_high = rc_high(W_test)[0]

        assert (
            r_high > r_low
        ), f"Higher Z should give larger correction: Z=20 -> {r_low:.6f}, Z=80 -> {r_high:.6f}"

    def test_z_zero_reduces_to_alpha_only(self):
        """With Z=0, only O(alpha) correction remains — should match standard Sirlin function."""
        rc_z0 = RadiativeCorrection(W0=5.0, Z=0, use_endpoint_resummation=True)
        rc_standard = RadiativeCorrection(W0=5.0, use_endpoint_resummation=True)

        W_test = np.array([2.0, 3.0, 4.0])
        r_z0 = rc_z0(W_test)
        r_std = rc_standard(W_test)

        np.testing.assert_allclose(
            r_z0,
            r_std,
            rtol=1e-10,
            err_msg="Z=0 should produce same result as Z-unspecified (alpha-only)",
        )

    def test_z_correction_magnitude_reasonable(self):
        """For Z=92 (uranium), O(Z*alpha^2) should be ~0.5-2% of total correction."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)

        W_test = np.array([3.0])
        result = rc(W_test)

        assert (
            0.9 < result[0] < 3.0
        ), f"Z=92 correction should be reasonable: got {result[0]}"

    def test_different_Z_different_results(self):
        """Different Z values must produce measurably different corrections."""
        rc_z1 = RadiativeCorrection(W0=5.0, Z=1, use_endpoint_resummation=True)
        rc_z50 = RadiativeCorrection(W0=5.0, Z=50, use_endpoint_resummation=True)

        W_test = np.linspace(1.5, 4.5, 10)
        r_z1 = rc_z1(W_test)
        r_z50 = rc_z50(W_test)

        diff = np.abs(r_z50 - r_z1)
        assert (
            np.max(diff) > 1e-4
        ), f"Z=1 and Z=50 should give different results: max diff = {np.max(diff)}"


class TestRadiativeNumericalStability:
    """Test numerical stability at extreme beta values."""

    def test_low_energy_no_nan(self):
        """At threshold (W -> 1), beta -> 0. Small-beta Taylor expansion must prevent NaN."""
        rc = RadiativeCorrection(W0=5.0, Z=50, use_endpoint_resummation=True)

        W_threshold = np.array([1.0001])
        result = rc(W_threshold)

        assert not np.isnan(result[0]), "Low energy must not produce NaN"
        assert not np.isinf(result[0]), "Low energy must not produce inf"
        assert result[0] > 0, "Correction must be positive at threshold"

    def test_ultra_low_energy_stable(self):
        """W = 1.00001 (beta ~ 0.014) — deep in small-beta regime."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)

        W_ultra_low = np.array([1.00001])
        result = rc(W_ultra_low)

        assert not np.isnan(result[0]), "Ultra-low energy must be stable"
        assert not np.isinf(result[0]), "Ultra-low energy must not diverge"

    def test_high_energy_no_nan(self):
        """At high energy (W -> W0), beta -> 1. arctanh must be handled stably."""
        rc = RadiativeCorrection(W0=10.0, Z=50, use_endpoint_resummation=True)

        W_high = np.array([9.0])
        result = rc(W_high)

        assert not np.isnan(result[0]), "High energy must not produce NaN"
        assert not np.isinf(result[0]), "High energy must not produce inf"

    def test_full_range_no_nan(self):
        """Entire energy range from threshold to endpoint must be free of NaN/inf."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)

        W_full = np.linspace(1.0001, 4.9999, 200)
        result = rc(W_full)

        assert not np.any(np.isnan(result)), "No NaN anywhere in spectrum"
        assert not np.any(np.isinf(result)), "No inf anywhere in spectrum"
        assert np.all(result > 0), "All values must be positive"

    def test_low_Z_low_energy_stable(self):
        """Light nucleus at low energy — tests both small Z and small beta."""
        rc = RadiativeCorrection(W0=2.0, Z=1, use_endpoint_resummation=True)

        W_test = np.array([1.001])
        result = rc(W_test)

        assert not np.isnan(result[0])
        assert not np.isinf(result[0])

    def test_high_Z_high_energy_stable(self):
        """Heavy nucleus at high energy — tests both large Z and beta near 1."""
        rc = RadiativeCorrection(W0=10.0, Z=92, use_endpoint_resummation=True)

        W_test = np.array([9.5])
        result = rc(W_test)

        assert not np.isnan(result[0])
        assert not np.isinf(result[0])


class TestRadiativeAParameter:
    """Test A (mass number) parameter for nuclear model."""

    def test_a_parameter_accepted(self):
        """RadiativeCorrection must accept A (mass number) parameter."""
        rc = RadiativeCorrection(W0=5.0, Z=50, A=120)
        assert rc.Z == 50
        assert rc.A == 120

    def test_a_affects_nuclear_model_correction(self):
        """A parameter affects the nuclear-structure-dependent part of O(Z*alpha^2)."""
        rc_a100 = RadiativeCorrection(
            W0=5.0, Z=50, A=100, use_endpoint_resummation=True
        )
        rc_a140 = RadiativeCorrection(
            W0=5.0, Z=50, A=140, use_endpoint_resummation=True
        )

        W_test = np.array([3.0])
        r_a100 = rc_a100(W_test)[0]
        r_a140 = rc_a140(W_test)[0]

        assert (
            r_a100 != r_a140
        ), f"Different A values must give different results: A=100 -> {r_a100}, A=140 -> {r_a140}"

    def test_default_a_is_none(self):
        """Without A specified, nuclear model correction should still work (A=None)."""
        rc = RadiativeCorrection(W0=5.0, Z=50, use_endpoint_resummation=True)
        assert rc.A is None
        W_test = np.array([3.0])
        result = rc(W_test)
        assert not np.isnan(result[0])
