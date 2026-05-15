"""
Code quality / API tests for FiniteSizeL0 and ChargeDistributionU.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.finite_size import FiniteSizeL0, ChargeDistributionU


class TestFiniteSizeL0API:
    @pytest.mark.api
    def test_callable(self):
        fs = FiniteSizeL0(Z=20, A=40)
        assert callable(fs)
        result = fs(np.array([2.0]))
        assert isinstance(result, np.ndarray)

    @pytest.mark.api
    def test_repr(self):
        fs = FiniteSizeL0(Z=20, A=40)
        assert "FiniteSizeL0" in repr(fs)

class TestFiniteSizeL0Properties:
    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [2.0],
        np.linspace(1.05, 6.0, 20),
    ])
    def test_positive(self, W):
        fs = FiniteSizeL0(Z=92, A=238)
        from tests.common.test_property_tests import _test_positive_values
        _test_positive_values(fs, W, component_name="FiniteSizeL0")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape(self, n):
        fs = FiniteSizeL0(Z=20, A=40)
        W = np.linspace(1.1, 3.0, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(fs, W, component_name="FiniteSizeL0")


class TestChargeDistributionUAPI:
    @pytest.mark.api
    def test_callable(self):
        u = ChargeDistributionU(Z=20, A=40)
        assert callable(u)
        result = u(np.array([2.0]))
        assert isinstance(result, np.ndarray)

    @pytest.mark.api
    def test_repr(self):
        u = ChargeDistributionU(Z=20, A=40)
        assert "ChargeDistributionU" in repr(u)


class TestChargeDistributionUProperties:
    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [2.0],
        np.linspace(1.1, 5.0, 10),
    ])
    def test_non_negative(self, W):
        u = ChargeDistributionU(Z=20, A=40)
        from tests.common.test_property_tests import _test_non_negative_values
        _test_non_negative_values(u, W, component_name="ChargeDistributionU")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape(self, n):
        u = ChargeDistributionU(Z=20, A=40)
        W = np.linspace(1.1, 3.0, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(u, W, component_name="ChargeDistributionU")
