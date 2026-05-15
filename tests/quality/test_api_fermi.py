"""
Code quality / API tests for FermiFunction.

Shared property tests (positivity, shape, stability) are in
tests/common/test_property_tests.py and are NOT duplicated here.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction


class TestFermiAPI:
    """Tests for the FermiFunction API contract."""

    @pytest.mark.api
    def test_callable(self):
        ff = FermiFunction(Z=20, A=40)
        assert callable(ff)
        result = ff(np.array([1.5]))
        assert isinstance(result, np.ndarray)

    @pytest.mark.api
    def test_repr(self):
        ff = FermiFunction(Z=20, A=40)
        assert "FermiFunction" in repr(ff)

class TestFermiProperties:
    """Shared property tests for FermiFunction."""

    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [1.5],
        [1.1, 2.0, 3.0],
        np.linspace(1.01, 5.0, 20),
    ])
    def test_positive(self, W):
        """F₀ must be positive everywhere — it's a probability enhancement factor."""
        ff = FermiFunction(Z=92, A=238)  # Uranium — extreme case
        from tests.common.test_property_tests import _test_positive_values
        _test_positive_values(ff, W, component_name="FermiFunction")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape_preserved(self, n):
        """Output shape must match input shape."""
        ff = FermiFunction(Z=20, A=40)
        W = np.linspace(1.1, 3.0, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(ff, W, component_name="FermiFunction")

    @pytest.mark.quality
    @pytest.mark.parametrize("W", [
        np.linspace(1.01, 6.0, 50),
        np.linspace(1.001, 10.0, 100),
    ])
    def test_no_nan(self, W):
        """For very heavy nuclei, the loggamma implementation must not overflow."""
        ff = FermiFunction(Z=92, A=238)
        from tests.common.test_property_tests import _test_no_nan_across_range
        _test_no_nan_across_range(ff, W, component_name="FermiFunction")

    @pytest.mark.quality
    def test_no_nan_at_threshold(self):
        """At W→1, p→0 and η→∞. The implementation must handle this limit."""
        ff = FermiFunction(Z=50, A=120)  # Tin — moderate heavy
        result = ff(np.array([1.001]))
        assert not np.isnan(result[0]), "Must handle threshold (W→1) without NaN"
        assert result[0] > 0, "Threshold value must be positive and finite"
