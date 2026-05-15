"""
Physics tests for ScreeningCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.screening import ScreeningCorrection


class TestScreeningPhysics:
    """Test expected physical behavior of the screening correction."""

    @pytest.mark.physics
    def test_vanishes_at_high_energy(self):
        """At high energies (W >> V₀), S(W) → 1. The logistic switch should be ~0 here."""
        ff = FermiFunction(Z=19, A=40)
        sc = ScreeningCorrection(ff)

        high_W = np.array([3.0, 5.0])
        result = sc(high_W)

        assert (
            abs(result[0] - 1.0) < 0.1
        ), f"Screening must approach unity at high energy: got {result}"

    @pytest.mark.physics
    def test_reduces_at_low_energy(self):
        """At low energy, screening reduces the effective Coulomb field → S(W) < 1.

        Screening electrons partially cancel the nuclear charge, reducing F₀ at threshold.
        """
        ff = FermiFunction(Z=30, A=70)
        sc = ScreeningCorrection(ff)

        low_W = np.array([1.05])
        mid_W = np.array([2.0])

        s_low = sc(low_W)[0]
        s_mid = sc(mid_W)[0]

        assert s_low < 1.1, f"Low-energy screening should suppress: {s_low}"
        assert abs(s_mid - 1.0) < abs(
            s_low - 1.0
        ), "Correction must decrease toward unity as W increases"

    @pytest.mark.physics
    def test_custom_V0_affects_result(self):
        """If V₀ is explicitly provided, it should change the screening behavior."""
        ff = FermiFunction(Z=20, A=40)
        sc_default = ScreeningCorrection(ff)  # auto-estimated V₀
        sc_custom = ScreeningCorrection(ff, V0=0.1)

        low_W = np.array([1.05])
        result_default = sc_default(low_W)[0]
        result_custom = sc_custom(low_W)[0]

        assert (
            abs(result_default - result_custom) > 0
        ), "Custom V₀ must change the screening behavior"
