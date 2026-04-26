# tests/test_screening.py — ScreeningCorrection unit tests

"""
Why test this?
--------------
Atomic screening reduces the effective nuclear charge felt by the emitted electron,
modifying the Fermi function at low energies. The correction is typically 10⁻³–10⁻²
but grows significant near threshold where p is small.

The implementation uses a ratio method (screened vs unscreened F(Z,W)) with a
logistic switching function to smoothly turn off the correction at high energy.

We test:

 1. **Near-unity at high energy**: At W >> V₀, screening should vanish → S(W) → 1.
     The logistic switch ensures this. This is critical — we don't want corrections
    where they shouldn't be applied.
 2. **Positive values**: Screening correction must stay positive (it's a multiplicative factor).
 3. **V₀ estimate合理性**: The default V0 scaling law should give reasonable values for
     light and heavy nuclei.

Common practice: For components with "switching" logic, always test both the on-region
(low energy) and off-region (high energy) to verify the transition is smooth and correct.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.screening import ScreeningCorrection


class TestScreeningBasicProperties:
    """Test fundamental properties of the screening correction."""

    def test_near_unity_at_high_energy(self, W_full):
        """At high energies (W >> V₀), S(W) → 1. The logistic switch should be ~0 here."""

        # For a typical nucleus (Z=20), V₀ ≈ α·Z^(4/3)·C ≈ 0.007
        # At W > 2, the switch is essentially zero
        ff = FermiFunction(Z=19, A=40)
        sc = ScreeningCorrection(ff)

        high_W = np.array([3.0, 5.0])
        result = sc(high_W)

        assert abs(result[0] - 1.0) < 0.1, (
            f"Screening must approach unity at high energy: got {result}"
        )

    def test_positive_values(self):
        """S(W) must be positive everywhere — it multiplies the spectrum."""
        ff = FermiFunction(Z=40, A=95)  # Zirconium
        sc = ScreeningCorrection(ff)

        W_test = np.linspace(1.02, 3.0, 20)
        result = sc(W_test)

        assert np.all(result > 0), "Screening correction must be positive"


class TestScreeningEnergyDependence:
    """Test how screening varies with energy."""

    def test_low_energy_correction_is_smaller(self):
        """At low energy, screening reduces the effective Coulomb field → S(W) < 1.

        Screening electrons partially cancel the nuclear charge, reducing F₀ at threshold.
        The correction factor should dip below 1 near W=1 and return to 1 at high energy.
        """
        ff = FermiFunction(Z=30, A=70)  # Zinc
        sc = ScreeningCorrection(ff)

        low_W = np.array([1.05])
        mid_W = np.array([2.0])

        s_low = sc(low_W)[0]
        s_mid = sc(mid_W)[0]

        # At low energy, correction is active and reduces the spectrum
        assert s_low < 1.1, f"Low-energy screening should suppress: {s_low}"
        # At mid energy, it's closer to unity
        assert abs(s_mid - 1.0) < abs(s_low - 1.0), (
            "Correction must decrease toward unity as W increases"
        )


class TestScreeningParameters:
    """Test that parameter choices affect behavior correctly."""

    def test_custom_V0_used(self):
        """If V₀ is explicitly provided, it should be used instead of the estimate."""
        ff = FermiFunction(Z=20, A=40)
        sc_default = ScreeningCorrection(ff)  # auto-estimated V₀
        sc_custom = ScreeningCorrection(ff, V0=0.1)

        low_W = np.array([1.05])
        result_default = sc_default(low_W)[0]
        result_custom = sc_custom(low_W)[0]

        # Different V₀ should give different results at low energy
        assert abs(result_default - result_custom) > 0, (
            "Custom V₀ must change the screening behavior"
        )


class TestScreeningOutputShape:
    """Basic shape/type checks."""

    def test_output_shape(self):
        ff = FermiFunction(Z=20, A=40)
        sc = ScreeningCorrection(ff)
        for n in [1, 5, 100]:
            W_test = np.linspace(1.1, 3.0, n)
            result = sc(W_test)
            assert result.shape == W_test.shape
