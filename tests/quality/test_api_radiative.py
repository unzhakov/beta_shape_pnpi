"""
Code quality / API tests for RadiativeCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.radiative import RadiativeCorrection


class TestRadiativeAPI:
    @pytest.mark.api
    def test_callable(self):
        rc = RadiativeCorrection(W0=5.0)
        assert callable(rc)
        result = rc(np.array([3.0]))
        assert isinstance(result, np.ndarray)

    @pytest.mark.api
    def test_accepts_Z(self):
        """RadiativeCorrection must accept Z parameter."""
        for Z in [1, 20, 50, 92]:
            rc = RadiativeCorrection(W0=5.0, Z=Z)
            assert rc.Z == Z

    @pytest.mark.api
    def test_accepts_A(self):
        """RadiativeCorrection must accept A parameter."""
        rc = RadiativeCorrection(W0=5.0, Z=50, A=120)
        assert rc.A == 120

    @pytest.mark.api
    def test_default_A_is_none(self):
        rc = RadiativeCorrection(W0=5.0, Z=50)
        assert rc.A is None


class TestRadiativeProperties:
    @pytest.mark.property
    @pytest.mark.parametrize("W", [
        [2.0],
        np.linspace(1.1, 4.9, 20),
    ])
    def test_positive(self, W):
        rc = RadiativeCorrection(W0=5.0)
        from tests.common.test_property_tests import _test_positive_values
        _test_positive_values(rc, W, component_name="RadiativeCorrection")

    @pytest.mark.quality
    @pytest.mark.parametrize("n", [1, 5, 100])
    def test_shape(self, n):
        rc = RadiativeCorrection(W0=5.0)
        W = np.linspace(1.1, 4.9, n)
        from tests.common.test_property_tests import _test_output_shape_preserved
        _test_output_shape_preserved(rc, W, component_name="RadiativeCorrection")

    @pytest.mark.quality
    @pytest.mark.parametrize("W0", [2.0, 5.0, 10.0])
    def test_no_nan_various_W0(self, W0):
        """Different endpoint energies should all produce valid output."""
        rc = RadiativeCorrection(W0=W0)
        W = np.linspace(1.1, min(9.9, W0 - 0.01), 20)
        from tests.common.test_property_tests import _test_no_nan_across_range
        _test_no_nan_across_range(rc, W, component_name=f"RadiativeCorrection(W0={W0})")

    @pytest.mark.quality
    @pytest.mark.parametrize("W", [
        np.linspace(1.0001, 4.9999, 200),
        np.linspace(1.001, 9.99, 100),
    ])
    def _test_full_range_stability(self, W):
        """Entire energy range from threshold to endpoint must be free of NaN/inf."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)
        from tests.common.test_property_tests import _test_full_range_stability
        _test_full_range_stability(rc, W, component_name="RadiativeCorrection(Z=92)")

    @pytest.mark.quality
    def test_ultra_low_energy_stable(self):
        """W = 1.00001 (beta ~ 0.014) — deep in small-beta regime."""
        rc = RadiativeCorrection(W0=5.0, Z=92, use_endpoint_resummation=True)
        W_ultra_low = np.array([1.00001])
        result = rc(W_ultra_low)
        assert not np.isnan(result[0])
        assert not np.isinf(result[0])
