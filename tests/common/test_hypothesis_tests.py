"""
Property-based tests using Hypothesis for edge case discovery.

While the parametrized property tests in test_property_tests.py cover
common cases, Hypothesis generates thousands of random inputs to find
edge cases that human-written tests might miss.

These tests are marked @pytest.mark.slow because they can take several
seconds to complete. They are intended for CI/full test runs, not for
fast feedback during development.
"""

from typing import List

import numpy as np
import pytest
from hypothesis import given, settings, strategies as st

from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.phase_space import PhaseSpace


# ---------------------------------------------------------------------------
# Fermi function: random Z, A, W combinations
# ---------------------------------------------------------------------------

@given(
    Z=st.integers(min_value=1, max_value=92),
    A=st.integers(min_value=1, max_value=300),
    W=st.floats(min_value=1.01, max_value=10.0, allow_infinity=False, allow_nan=False),
)
@settings(max_examples=200, deadline=None)
def test_fermi_stable_for_random_params(Z: int, A: int, W: float):
    """Fermi function must not crash or produce NaN for any valid (Z, A, W)."""
    ff = FermiFunction(Z=Z, A=A)
    result = ff(np.array([W]))

    assert not np.isnan(result[0]), f"NaN for Z={Z}, A={A}, W={W}"
    assert not np.isinf(result[0]), f"Inf for Z={Z}, A={A}, W={W}"
    assert result[0] > 0, f"Non-positive for Z={Z}, A={A}, W={W}: {result[0]}"


@given(
    Z=st.integers(min_value=10, max_value=92),
)
@settings(max_examples=200, deadline=None)
def test_fermi_monotonic_for_fixed_Z(Z: int):
    """For high Z, Fermi function should be monotonically decreasing with W.

    For Z=1 the finite-size correction can cause slight non-monotonicity,
    so we test high-Z cases where Coulomb dominates.
    """
    ff = FermiFunction(Z=Z, A=Z * 2 + 10)
    W_arr = np.linspace(1.05, 4.0, 10)
    result = ff(W_arr)

    diffs = np.diff(result)
    assert np.all(diffs <= 0), (
        f"Fermi function not monotonically decreasing for Z={Z}: "
        f"diffs = {diffs}"
    )


# ---------------------------------------------------------------------------
# Phase space: random W0, W combinations
# ---------------------------------------------------------------------------

@given(
    W0=st.floats(min_value=2.0, max_value=10.0, allow_infinity=False, allow_nan=False),
    W=st.floats(min_value=1.01, max_value=5.0, allow_infinity=False, allow_nan=False),
)
@settings(max_examples=200, deadline=None)
def test_phase_space_valid_for_random_params(W0: float, W: float):
    """Phase space must not crash for any valid (W0, W) with W < W0."""
    # Ensure W is always less than W0
    W = min(W, W0 - 0.02)
    ps = PhaseSpace(W0=W0)
    result = ps(np.array([W]))

    assert not np.isnan(result[0]), f"NaN for W0={W0:.2f}, W={W:.2f}"
    # Phase space should be non-negative for W < W0
    assert result[0] >= 0, f"Negative for W0={W0:.2f}, W={W:.2f}: {result[0]}"


@given(
    W0=st.floats(min_value=2.0, max_value=10.0, allow_infinity=False, allow_nan=False),
)
@settings(max_examples=100, deadline=None)
def test_phase_space_threshold_at_W1(W0: float):
    """Phase space must be exactly zero at W=1."""
    ps = PhaseSpace(W0=W0)
    result = ps(np.array([1.0]))

    assert result[0] == 0.0, f"Phase space should be 0 at W=1 for W0={W0:.2f}, got {result[0]}"


@given(
    W0=st.floats(min_value=2.0, max_value=10.0, allow_infinity=False, allow_nan=False),
)
@settings(max_examples=100, deadline=None)
def test_phase_space_endpoint_at_W0(W0: float):
    """Phase space must be exactly zero at W=W0."""
    ps = PhaseSpace(W0=W0)
    result = ps(np.array([W0]))

    assert result[0] == 0.0, f"Phase space should be 0 at W=W0={W0:.2f}, got {result[0]}"
