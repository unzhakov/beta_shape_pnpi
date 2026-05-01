# tests/test_finite_size.py — FiniteSizeL0 and ChargeDistributionU unit tests

"""
Why test this?
--------------
Finite nuclear size corrections are small (typically 10⁻³–10⁻²) but systematic.
They depend on the nuclear radius R ∝ A^(1/3), so they must scale correctly
with atomic mass number A.

We verify:

 1. **L₀ starts near unity**: For low-Z nuclei, L₀ ≈ 1 + small correction.
     The correction is O(αZ·WR) which is tiny for light elements.
 2. **Charge distribution U term > 0**: This is a second-order positive term
     proportional to (αZW R)² — always enhances the spectrum.
 3. **A-scaling**: Larger A → larger radius → larger correction.
 4. **Combined effect ≈ unity for light nuclei**: For Z=6, A=12, the total
     finite-size correction should be very close to 1 (well within 5%).

Common practice: For small-correction components like this, test that they are
"close enough to 1" for boundary cases where the physics says the effect is negligible.
This catches sign errors and magnitude bugs early.
"""

import numpy as np

from beta_spectrum.components.finite_size import FiniteSizeL0, ChargeDistributionU


class TestFiniteSizeL0:
    """Tests for the L₀ finite nuclear size correction."""

    def test_near_unity_for_light_nucleus(self):
        """For Z=6 (carbon), L₀ should be very close to 1 — within a few percent.

        The correction is O(αZ·WR) which for carbon is ~0.001. If we get >>5%,
         something is wrong with the expansion coefficients or sign convention.
        """
        fs = FiniteSizeL0(Z=6, A=12)
        W_test = np.array([2.0])  # Mid-range energy
        result = fs(W_test)

        assert (
            abs(result[0] - 1.0) < 0.05
        ), f"L₀ for Z=6 should be near unity but got {result[0]}"

    def test_correction_increases_with_Z(self):
        """Higher Z → larger αZ term in L₀ expansion → larger deviation from 1."""
        fs_z6 = FiniteSizeL0(Z=6, A=12)
        fs_z30 = FiniteSizeL0(Z=30, A=70)
        fs_z80 = FiniteSizeL0(Z=80, A=200)

        W_test = np.array([2.5])

        dev_6 = abs(fs_z6(W_test)[0] - 1.0)
        dev_30 = abs(fs_z30(W_test)[0] - 1.0)
        dev_80 = abs(fs_z80(W_test)[0] - 1.0)

        assert (
            dev_80 > dev_30 > dev_6
        ), f"Deviation from unity must increase with Z: {dev_6}, {dev_30}, {dev_80}"

    def test_positive_everywhere(self):
        """L₀ should be positive — it's a correction factor multiplicative to the spectrum."""
        fs = FiniteSizeL0(Z=92, A=238)
        W_test = np.linspace(1.05, 6.0, 20)
        result = fs(W_test)

        assert np.all(result > 0), "Finite size L₀ must be positive"


class TestChargeDistributionU:
    """Tests for the U charge distribution correction."""

    def test_near_unity_for_light_nucleus(self):
        """For Z=6, A term is tiny — (αZ·WR)²/5 ≈ 10⁻⁶. Should be essentially 1."""
        u = ChargeDistributionU(Z=6, A=12)
        result = u(np.array([2.0]))

        assert (
            abs(result[0] - 1.0) < 0.01
        ), f"U for Z=6 should be near unity but got {result[0]}"

    def test_always_positive_correction(self):
        """U = 1 + (1/5)(αZW R)² is always > 1 since the squared term is positive."""
        u = ChargeDistributionU(Z=20, A=40)
        W_test = np.linspace(1.1, 5.0, 10)
        result = u(W_test)

        assert np.all(result >= 1.0), "U must be ≥ 1 for all W (squared term)"

    def test_increases_with_Z(self):
        """Higher Z → larger U correction."""
        u_z20 = ChargeDistributionU(Z=20, A=40)
        u_z60 = ChargeDistributionU(Z=60, A=150)

        W_test = np.array([3.0])

        assert (
            u_z60(W_test)[0] > u_z20(W_test)[0]
        ), f"U must increase with Z: {u_z20(W_test)} < {u_z60(W_test)}"


class TestFiniteSizeCombined:
    """Test the combined effect of L₀ and U corrections."""

    def test_combined_near_unity_for_carbon(self):
        """For carbon-12, both corrections are tiny. The product should be ~1."""
        l0 = FiniteSizeL0(Z=6, A=12)
        u = ChargeDistributionU(Z=6, A=12)

        W_test = np.linspace(1.5, 4.0, 10)
        combined = l0(W_test) * u(W_test)

        assert np.all(
            np.abs(combined - 1.0) < 0.05
        ), f"Combined finite-size for carbon should be near unity: {combined}"

    def test_combined_positive(self):
        """Product of two positive factors must be positive."""
        l0 = FiniteSizeL0(Z=92, A=238)
        u = ChargeDistributionU(Z=92, A=238)

        W_test = np.linspace(1.05, 6.0, 20)
        combined = l0(W_test) * u(W_test)

        assert np.all(combined > 0), "Combined correction must be positive"


class TestFiniteSizeOutputShape:
    """Basic shape checks."""

    def test_output_shape(self):
        fs = FiniteSizeL0(Z=20, A=40)
        for n in [1, 5, 100]:
            W_test = np.linspace(1.1, 3.0, n)
            result = fs(W_test)
            assert result.shape == W_test.shape

    def test_u_output_shape(self):
        u = ChargeDistributionU(Z=20, A=40)
        for n in [1, 5, 100]:
            W_test = np.linspace(1.1, 3.0, n)
            result = u(W_test)
            assert result.shape == W_test.shape
