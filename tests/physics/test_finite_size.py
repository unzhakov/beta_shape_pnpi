"""
Physics tests for FiniteSizeL0 and ChargeDistributionU.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.finite_size import FiniteSizeL0, ChargeDistributionU


class TestFiniteSizeL0Physics:
    """Tests for the L₀ finite nuclear size correction."""

    @pytest.mark.physics
    def test_near_unity_for_light_nucleus(self):
        """For Z=6 (carbon), L₀ should be very close to 1 — within a few percent."""
        fs = FiniteSizeL0(Z=6, A=12)
        W_test = np.array([2.0])
        result = fs(W_test)

        assert (
            abs(result[0] - 1.0) < 0.05
        ), f"L₀ for Z=6 should be near unity but got {result[0]}"

    @pytest.mark.physics
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


class TestChargeDistributionUPhysics:
    """Tests for the U charge distribution correction."""

    @pytest.mark.physics
    def test_near_unity_for_light_nucleus(self):
        """For Z=6, U is essentially 1."""
        u = ChargeDistributionU(Z=6, A=12)
        result = u(np.array([2.0]))
        assert abs(result[0] - 1.0) < 0.01

    @pytest.mark.physics
    def test_always_above_unity(self):
        """U = 1 + (1/5)(αZW R)² is always ≥ 1 since the squared term is positive."""
        u = ChargeDistributionU(Z=20, A=40)
        W_test = np.linspace(1.1, 5.0, 10)
        result = u(W_test)
        assert np.all(result >= 1.0), "U must be ≥ 1 for all W (squared term)"

    @pytest.mark.physics
    def test_increases_with_Z(self):
        """Higher Z → larger U correction."""
        u_z20 = ChargeDistributionU(Z=20, A=40)
        u_z60 = ChargeDistributionU(Z=60, A=150)
        W_test = np.array([3.0])
        assert u_z60(W_test)[0] > u_z20(W_test)[0]


class TestFiniteSizeCombinedPhysics:
    """Combined L₀ × U physics tests."""

    @pytest.mark.physics
    def test_combined_near_unity_for_carbon(self):
        """For carbon-12, both corrections are tiny. The product should be ~1."""
        l0 = FiniteSizeL0(Z=6, A=12)
        u = ChargeDistributionU(Z=6, A=12)

        W_test = np.linspace(1.5, 4.0, 10)
        combined = l0(W_test) * u(W_test)

        assert np.all(
            np.abs(combined - 1.0) < 0.05
        ), f"Combined finite-size for carbon should be near unity: {combined}"
