"""
Physics tests for FermiFunction — unique physical behavior.

Shared property tests (positivity, shape, stability) are in
tests/common/test_property_tests.py and are NOT duplicated here.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction


class TestFermiPhysicalBehavior:
    """Test expected physical behavior of the Fermi function."""

    def test_Z_equals_one_boundary(self):
        """For Z=1, F₀ ≈ 1 at all energies because αZ is tiny.

        This is the most important sanity check. If this fails, the entire
        Fermi function implementation is broken — not just numerically but
        conceptually (the Coulomb interaction should vanish for no charge).
        """
        ff = FermiFunction(Z=1, A=1)
        W_test = np.array([1.5, 2.0, 3.0, 5.0])
        result = ff(W_test)

        # For Z=1, the correction from finite nuclear size and Coulomb terms gives ~3% deviation
        assert np.all(
            np.abs(result - 1.0) < 0.05
        ), f"Z=1 Fermi function should be ≈1 but got {result}"

    @pytest.mark.physics
    def test_decreases_with_energy(self):
        """For β⁻ decay (positive Z_daughter), F₀ decreases as W increases.

        Reason: η = αZW/p. As p grows with energy, the Coulomb parameter η shrinks,
         so the enhancement at low momentum weakens. This is universal for all Z > 0.
        """
        ff = FermiFunction(Z=20, A=40)
        W_test = np.linspace(1.1, 2.0, 5)
        result = ff(W_test)

        # Check monotonic decrease (each successive point should be smaller)
        diffs = np.diff(result)
        assert np.all(
            diffs < 0
        ), f"Fermi function must decrease with W; got diffs {diffs}"

    @pytest.mark.physics
    def test_Z_scaling(self):
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

        assert (
            f_80 > f_30 > f_10
        ), f"Fermi function must increase with Z: {f_10} < {f_30} < {f_80}"

    @pytest.mark.physics
    def test_uranium_threshold_enhancement(self):
        """For uranium (Z=92), F₀ at threshold should be orders of magnitude
        larger than at mid-energy due to strong Coulomb attraction.
        """
        ff = FermiFunction(Z=92, A=238)

        W_low = np.array([1.01])   # Near threshold
        W_mid = np.array([3.0])     # Mid-range

        f_low = ff(W_low)[0]
        f_mid = ff(W_mid)[0]

        assert f_low > f_mid * 10, (
            f"Uranium threshold enhancement should be large: "
            f"F(1.01)={f_low:.1f} vs F(3.0)={f_mid:.1f}"
        )
