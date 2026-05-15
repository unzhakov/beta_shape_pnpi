"""
Code quality / API tests for ExchangeCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.exchange import ExchangeCorrection


class TestExchangeAPI:
    @pytest.mark.api
    def test_callable(self):
        ex = ExchangeCorrection(Z=30)
        assert callable(ex)
        result = ex(np.array([2.0]))
        assert isinstance(result, np.ndarray)

class TestExchangeProperties:
    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [2.0],
        np.linspace(1.01, 6.0, 50),
    ])
    def test_positive(self, W):
        ex = ExchangeCorrection(Z=92)
        from tests.common.test_property_tests import _test_positive_values
        _test_positive_values(ex, W, component_name="ExchangeCorrection")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape(self, n):
        ex = ExchangeCorrection(Z=30)
        W = np.linspace(1.1, 3.0, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(ex, W, component_name="ExchangeCorrection")

    @pytest.mark.quality
    @pytest.mark.parametrize("W", [
        np.linspace(1.01, 6.0, 50),
    ])
    def test_no_nan(self, W):
        ex = ExchangeCorrection(Z=92)
        from tests.common.test_property_tests import _test_no_nan_across_range
        _test_no_nan_across_range(ex, W, component_name="ExchangeCorrection")
