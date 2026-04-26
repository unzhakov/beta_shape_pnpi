# tests/test_fermi.py — FermiFunction component unit tests

"""
Why test this?
--------------
The Fermi function F₀(Z,W) is the most critical correction for heavy nuclei.
At low Z it's ~1, but at high Z (e.g., Z=92 uranium) it can reach 10⁴+ near
threshold because the Coulomb attraction pulls electrons toward the nucleus.

We test:

 1. **Z=1 sanity**: For hydrogen, F₀ ≈ 1 (no Coulomb distortion). This is a
    well-known boundary condition — if Z=1 doesn't give ~1, something fundamental
    is broken in the loggamma computation.
 2. **Monotonicity with energy**: For β⁻ decay (Z_daughter > 0), F₀ decreases
     with increasing W because p grows and η = αZW/p shrinks. This is a physical
    expectation: Coulomb effects are strongest at low momentum.
 3. **Z-scaling**: Higher Z → larger correction. F₀(Z=50) should >> F₀(Z=1).
 4. **Numerical stability**: The loggamma-based implementation avoids overflow
     that would occur with naive exp() of large arguments. We test this at high Z.

Common practice: When testing a function with special mathematical functions
(loggamma, spence, etc.), always include a "sanity check" against a known case
(Z=1) to catch implementation bugs early.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction


class TestFermiZeros:
    """Test Z=1 (hydrogen) boundary condition."""

    def test_Z_equals_one_gives_near_unity(self):
        """For Z=1, F₀ ≈ 1 at all energies because αZ is tiny.

        This is the most important sanity check. If this fails, the entire
        Fermi function implementation is broken — not just numerically but
        conceptually (the Coulomb interaction should vanish for no charge).
        """
        ff = FermiFunction(Z=1, A=1)
        W_test = np.array([1.5, 2.0, 3.0, 5.0])
        result = ff(W_test)

        # For Z=1, the correction from finite nuclear size and Coulomb terms gives ~3% deviation
        assert np.all(np.abs(result - 1.0) < 0.05), (
            f"Z=1 Fermi function should be ≈1 but got {result}"
        )


class TestFermiPhysicalBehavior:
    """Test expected physical behavior of the Fermi function."""

    def test_decreases_with_energy(self, W_mid):
        """For β⁻ decay (positive Z_daughter), F₀ decreases as W increases.

        Reason: η = αZW/p. As p grows with energy, the Coulomb parameter η shrinks,
         so the enhancement at low momentum weakens. This is universal for all Z > 0.
        """
        ff = FermiFunction(Z=20, A=40)
        result = ff(W_mid)

        # Check monotonic decrease (each successive point should be smaller)
        diffs = np.diff(result)
        assert np.all(diffs < 0), f"Fermi function must decrease with W; got diffs {diffs}"

    def test_increases_with_Z(self):
        """Higher Z → stronger Coulomb attraction → larger F₀ at low energy.

        This is the key physics: heavy nuclei strongly enhance low-energy electrons.
        We compare at a fixed mid-range energy (W=2) across three Z values.
        """
        ff_z10 = FermiFunction(Z=10, A=25)
        ff_z30 = FermiFunction(Z=30, A=65)
        ff_z80 = FermiFunction(Z=80, A=200)

        W_test = np.array([2.0])

        f_10 = ff_z10(W_test)[0]
        f_30 = ff_z30(W_test)[0]
        f_80 = ff_z80(W_test)[0]

        assert f_80 > f_30 > f_10, (
            f"Fermi function must increase with Z: {f_10} < {f_30} < {f_80}"
        )

    def test_positive_values(self):
        """F₀ must be positive everywhere — it's a probability enhancement factor."""
        ff = FermiFunction(Z=92, A=238)  # Uranium — extreme case
        W_test = np.linspace(1.01, 5.0, 20)
        result = ff(W_test)

        assert np.all(result > 0), "Fermi function must be positive for all Z and W"


class TestFermiNumericalStability:
    """Test that the loggamma implementation doesn't overflow or produce NaN."""

    def test_no_nan_at_high_Z(self):
        """For very heavy nuclei, naive computation of |Γ(γ+iη)|² would overflow.
        The loggamma trick is essential here — we verify it produces valid output.
        """
        ff = FermiFunction(Z=92, A=238)  # Uranium-238
        W_test = np.linspace(1.05, 6.0, 50)
        result = ff(W_test)

        assert not np.any(np.isnan(result)), "No NaN values allowed at high Z"
        assert not np.any(np.isinf(result)), "No infinite values allowed at high Z"

    def test_no_nan_at_threshold(self):
        """At W→1, p→0 and η→∞. The implementation must handle this limit."""
        ff = FermiFunction(Z=50, A=120)  # Tin — moderate heavy
        result = ff(np.array([1.001]))

        assert not np.isnan(result[0]), "Must handle threshold (W→1) without NaN"
        assert result[0] > 0, "Threshold value must be positive and finite"


class TestFermiOutputShape:
    """Test that output shape matches input — a basic correctness check."""

    def test_output_shape_preserved(self):
        ff = FermiFunction(Z=20, A=40)
        for n_points in [1, 5, 100]:
            W_test = np.linspace(1.1, 3.0, n_points)
            result = ff(W_test)
            assert result.shape == W_test.shape, (
                f"Shape mismatch: input {W_test.shape} → output {result.shape}"
            )
