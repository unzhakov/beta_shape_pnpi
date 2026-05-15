"""
Code quality / API tests for ScreeningCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.screening import ScreeningCorrection


class TestScreeningAPI:
    @pytest.mark.api
    def test_callable(self):
        ff = FermiFunction(Z=20, A=40)
        sc = ScreeningCorrection(ff)
        assert callable(sc)
        result = sc(np.array([2.0]))
        assert isinstance(result, np.ndarray)

    @pytest.mark.api
    def test_repr(self):
        ff = FermiFunction(Z=20, A=40)
        sc = ScreeningCorrection(ff)
        assert "Screening" in repr(sc)

class TestScreeningProperties:
    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [2.0],
        np.linspace(1.02, 3.0, 20),
    ])
    def test_positive(self, W):
        ff = FermiFunction(Z=40, A=95)
        sc = ScreeningCorrection(ff)
        from tests.common.test_property_tests import _test_positive_values
        _test_positive_values(sc, W, component_name="ScreeningCorrection")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape(self, n):
        ff = FermiFunction(Z=20, A=40)
        sc = ScreeningCorrection(ff)
        W = np.linspace(1.1, 3.0, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(sc, W, component_name="ScreeningCorrection")
