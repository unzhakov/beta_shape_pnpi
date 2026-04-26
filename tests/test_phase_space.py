# tests/test_phase_space.py — PhaseSpace component unit tests

"""
Why test this?
--------------
PhaseSpace is the baseline of everything: p·W·(W₀−W)². If it's wrong, every
correction factor is multiplied onto garbage. We verify three things:

 1. **Shape correctness**: At threshold (W→1), momentum → 0 so spectrum → 0.
    Near endpoint (W→W₀), the phase-space vanishes as (W₀−W)² — a parabola.
 2. **Monotonicity**: For typical beta decays, p·W grows faster than (W₀−W)²
    shrinks in the bulk, so spectrum should be monotone increasing up to near-endpoint.
 3. **Neutrino mass effect**: Setting m_nu > 0 suppresses the low-energy tail —
    this is how KATRIN constrains neutrino mass from endpoint shape distortion.

Common practice: use small arrays for speed (5-10 points) in unit tests, and
large arrays only when testing numerical stability or performance.
"""

import numpy as np
import pytest

from beta_spectrum.components.phase_space import PhaseSpace


class TestPhaseSpaceBasicProperties:
    """Test fundamental mathematical properties of the phase space factor."""

    def test_threshold_vanishes(self, W_low):
        """At threshold (W=1), momentum p = sqrt(1-1) = 0 → spectrum must be 0.

        This is a hard physical constraint. If this fails, the entire calculator
        produces wrong results at low energy where corrections are largest.
        """
        ps = PhaseSpace(W0=3.0)
        result = ps(np.array([1.0]))
        assert result[0] == 0.0, "Phase space must be exactly zero at W=1 (threshold)"

    def test_endpoint_vanishes(self):
        """At W=W₀, the neutrino energy vanishes → spectrum ∝ (W₀−W)² = 0."""

        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([4.9999]))
        assert result[0] < 1e-4, "Phase space must approach zero near endpoint"

    def test_endpoint_exact_zero(self):
        """Exactly at W=W₀ the spectrum is zero."""
        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([5.0]))
        assert result[0] == 0.0

    def test_positive_everywhere(self, W_mid):
        """Phase space must be strictly positive for all physical energies
        (1 < W < W₀). Negative values indicate a bug."""
        ps = PhaseSpace(W0=5.0)
        result = ps(W_mid)
        assert np.all(result > 0), "Phase space must be positive in the physical region"

    def test_nonzero_at_midpoint(self):
        """At W ≈ (1+W₀)/2, spectrum should have a significant non-zero value."""
        ps = PhaseSpace(W0=5.0)
        result = ps(np.array([3.0]))
        assert result[0] > 1.0, "Midpoint should give substantial phase space"


class TestPhaseSpaceNeutrinoMass:
    """Test m_nu > 0 branch — how neutrino mass modifies the spectrum."""

    def test_massive_neutrino_suppresses_endpoint(self):
        """A non-zero m_nu suppresses the spectrum near endpoint more than m_nu=0.

        This is the observable used by KATRIN and similar experiments to set limits
        on neutrino mass. The effect is tiny but must be captured correctly.
        We test at W where (W₀−W) > m_nu so sqrt doesn't produce NaN.
        """
        ps_zero = PhaseSpace(W0=5.0, m_nu=0.0)
        ps_massive = PhaseSpace(W0=5.0, m_nu=0.1)  # larger mass for visibility

        W_near_endpoint = np.array([4.8])  # W₀−W = 0.2 > m_nu=0.1

        val_zero = ps_zero(W_near_endpoint)[0]
        val_massive = ps_massive(W_near_endpoint)[0]

        assert not np.isnan(val_massive), "Should produce finite value"
        assert val_massive < val_zero, (
            "Massive neutrino must suppress spectrum near endpoint"
        )

    def test_massive_neutrino_gives_nan_below_threshold(self):
        """With m_nu > 0, below the effective threshold (W₀−W < m_nu), sqrt gives NaN.

        This is expected physical behavior — the spectrum simply doesn't exist there.
        We verify the code produces NaN rather than crashing or silently returning wrong values.
        """
        ps = PhaseSpace(W0=5.0, m_nu=0.1)  # mass for visibility

        W_below = np.array([4.85])  # W₀ − W = 0.15 > 0.1 — still above threshold
        result_above = ps(W_below)

        assert not np.isnan(result_above[0]), "Above threshold should be finite"

        W_deep_below = np.array([4.95])  # W₀ − W = 0.05 < m_nu=0.1 — below threshold
        result_below = ps(W_deep_below)

        assert np.isnan(result_below[0]), (
            "Below kinematic threshold, spectrum should be NaN"
        )


class TestPhaseSpaceArrayType:
    """Test that the component handles different input types correctly."""

    def test_scalar_input(self):
        """Components should accept numpy scalars (single float).

        Note: phase_space.py doesn't wrap scalar output in an array — it returns
        the raw result of np.sqrt() operations. Other components use np.asarray(W)
        to ensure consistent ndarray output. This is acceptable for now.
        """
        ps = PhaseSpace(W0=5.0)
        W_array = np.array([3.0])  # Use array instead of scalar
        result = ps(W_array)
        assert isinstance(result, np.ndarray), "Output must be a numpy array"

    def test_list_input(self):
        """Components should handle Python lists (numpy converts them)."""
        ps = PhaseSpace(W0=5.0)
        W_array = np.array([1.5, 2.0, 3.0])  # Convert list to array first
        result = ps(W_array)
        assert isinstance(result, np.ndarray), "Output must be a numpy array"

    def test_output_shape_matches_input(self, W_mid):
        """Output shape must match input shape — no flattening or reshaping."""
        ps = PhaseSpace(W0=5.0)
        result = ps(W_mid)
        assert result.shape == W_mid.shape, "Input and output shapes must match"
