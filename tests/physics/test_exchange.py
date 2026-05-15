"""
Physics tests for ExchangeCorrection.

Shared property tests in tests/common/test_property_tests.py.
"""

import numpy as np
import pytest

from beta_spectrum.components.exchange import ExchangeCorrection


class TestExchangePhysics:
    """Test expected physical behavior of the exchange correction."""

    @pytest.mark.physics
    def test_large_correction_heavy_nucleus(self):
        """For Z=80 (mercury), X >> 1 near threshold — can be > 20%.

        This is the defining characteristic of the exchange correction.
        """
        ex = ExchangeCorrection(Z=80)
        low_W = np.array([1.1])
        result = ex(low_W)
        assert result[0] > 1.0, f"Exchange must enhance spectrum at low energy for heavy Z: {result[0]}"

    @pytest.mark.physics
    def test_grows_as_energy_decreases(self):
        """X(W) should increase as W → 1 (lower momentum → stronger exchange)."""
        ex = ExchangeCorrection(Z=50)

        low_W = np.array([1.2])
        lower_W = np.array([1.05])

        s_low = ex(low_W)[0]
        s_lower = ex(lower_W)[0]

        assert s_lower > s_low, f"Exchange must increase at lower energy: {s_low} → {s_lower}"

    @pytest.mark.physics
    def test_vanishes_at_high_energy(self):
        """At W >> 1, the emitted electron is far from atomic electrons → X → 1."""
        ex = ExchangeCorrection(Z=40)
        high_W = np.array([3.0, 5.0])
        result = ex(high_W)
        assert abs(result[0] - 1.0) < 0.5, f"Exchange must approach unity at high energy: {result}"


class TestExchangeErrorHandling:
    """Test error handling for invalid inputs."""

    @pytest.mark.api
    def test_missing_Z_raises_error(self):
        """Z=1 has no electrons to exchange with — should fail gracefully."""
        with pytest.raises(ValueError, match="No exchange coefficients"):
            ExchangeCorrection(Z=1)
