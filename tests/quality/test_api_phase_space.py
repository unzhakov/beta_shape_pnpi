"""
Code quality / API tests for PhaseSpace.

Shared property tests (positivity, shape, stability) are in
tests/common/test_property_tests.py and are NOT duplicated here.
"""

import numpy as np
import pytest

from beta_spectrum.components.phase_space import PhaseSpace


class TestPhaseSpaceAPI:
    """Tests for the PhaseSpace API contract."""

    @pytest.mark.api
    def test_callable(self):
        ps = PhaseSpace(W0=5.0)
        assert callable(ps)
        result = ps(np.array([3.0]))
        assert isinstance(result, np.ndarray)

class TestPhaseSpaceProperties:
    """Shared property tests for PhaseSpace."""

    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [1.5],
        [1.1, 2.0, 3.0],
        np.linspace(1.01, 4.9, 20),
    ])
    def test_non_negative(self, W):
        """Phase space must be non-negative (zero at boundaries)."""
        ps = PhaseSpace(W0=5.0)
        from tests.common.test_property_tests import _test_non_negative_values
        _test_non_negative_values(ps, W, component_name="PhaseSpace")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape_preserved(self, n):
        """Output shape must match input shape."""
        ps = PhaseSpace(W0=5.0)
        W = np.linspace(1.1, 4.9, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(ps, W, component_name="PhaseSpace")

    @pytest.mark.quality
    @pytest.mark.parametrize("W", [
        np.linspace(1.01, 4.99, 50),
        np.linspace(1.001, 4.999, 100),
    ])
    def test_stability(self, W):
        """Phase space must be stable across the full energy range."""
        ps = PhaseSpace(W0=5.0)
        from tests.common.test_property_tests import _test_full_range_stability
        _test_full_range_stability(ps, W, component_name="PhaseSpace")
