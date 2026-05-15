"""
Physics tests for PhaseSpace — unique physical behavior.

Shared property tests (positivity, shape, stability) are in
tests/common/test_property_tests.py and are NOT duplicated here.
"""

import numpy as np
import pytest

from beta_spectrum.components.phase_space import PhaseSpace


class TestPhaseSpaceBasicProperties:
    """Test fundamental mathematical properties of the phase space factor."""

    @pytest.mark.physics
    def test_threshold_vanishes(self):
        """At threshold (W=1), momentum p = sqrt(1-1) = 0 → spectrum must be 0.

        This is a hard physical constraint. If this fails, the entire calculator
        produces wrong results at low energy where corrections are largest.
        """
        ps = PhaseSpace(W0=3.0)
        result = ps(np.array([1.0]))
        assert result[0] == 0.0, "Phase space must be exactly zero at W=1 (threshold)"

    @pytest.mark.physics
    def test_endpoint_vanishes(self):
        """At W=W₀, the neutrino energy vanishes → spectrum ∝ (W₀−W)² = 0."""
        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([4.9999]))
        assert result[0] < 1e-4, "Phase space must approach zero near endpoint"

    @pytest.mark.physics
    def test_endpoint_exact_zero(self):
        """Exactly at W=W₀ the spectrum is zero."""
        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([5.0]))
        assert result[0] == 0.0

    @pytest.mark.physics
    def test_positive_in_physical_range(self):
        """Phase space must be strictly positive for all physical energies (1 < W < W₀)."""
        ps = PhaseSpace(W0=5.0)
        W_test = np.linspace(1.1, 4.9, 10)
        result = ps(W_test)
        assert np.all(result > 0), "Phase space must be positive in the physical region"

    @pytest.mark.physics
    def test_midpoint_significant(self):
        """At W ≈ (1+W₀)/2, spectrum should have a significant non-zero value."""
        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([3.0]))
        assert result[0] > 1.0, "Midpoint should give substantial phase space"


class TestPhaseSpaceNeutrinoMass:
    """Test m_nu > 0 branch — how neutrino mass modifies the spectrum."""

    @pytest.mark.physics
    def test_massive_neutrino_suppresses_endpoint(self):
        """A non-zero m_nu suppresses the spectrum near endpoint more than m_nu=0.

        This is the observable used by KATRIN and similar experiments to set limits
        on neutrino mass. The effect is tiny but must be captured correctly.
        """
        ps_zero = PhaseSpace(W0=5.0, m_nu=0.0)
        ps_massive = PhaseSpace(W0=5.0, m_nu=0.1)  # larger mass for visibility

        W_near_endpoint = np.array([4.8])  # W₀−W = 0.2 > m_nu=0.1

        val_zero = ps_zero(W_near_endpoint)[0]
        val_massive = ps_massive(W_near_endpoint)[0]

        assert not np.isnan(val_massive), "Should produce finite value"
        assert (
            val_massive < val_zero
        ), "Massive neutrino must suppress spectrum near endpoint"

    @pytest.mark.physics
    def test_massive_neutrino_below_threshold(self):
        """With m_nu > 0, below the effective threshold (W₀−W < m_nu), sqrt gives NaN.

        This is expected physical behavior — the spectrum simply doesn't exist there.
        """
        ps = PhaseSpace(W0=5.0, m_nu=0.1)

        W_above = np.array([4.85])  # W₀ − W = 0.15 > 0.1 — still above threshold
        result_above = ps(W_above)

        assert not np.isnan(result_above[0]), "Above threshold should be finite"

        W_deep_below = np.array([4.95])  # W₀ − W = 0.05 < m_nu=0.1 — below threshold
        result_below = ps(W_deep_below)

        assert np.isnan(
            result_below[0]
        ), "Below kinematic threshold, spectrum should be NaN"
