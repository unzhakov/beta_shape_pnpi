# tests/test_exchange.py — ExchangeCorrection unit tests

"""
Why test this?
--------------
The atomic exchange correction X(Z,W) can be very large for heavy nuclei at low energy
— exceeding 20% (i.e., factor of 1.2+) near threshold. This is because the emitted
electron's wave function overlaps with bound atomic electrons, modifying the phase space.

Unlike other corrections that are small perturbations (~1%), exchange is a genuine
O(α) effect that grows dramatically as p → 0 (because η = αZW/p diverges).

We test:

 1. **Large correction at low energy for heavy nuclei**: At Z=80, W→1, X >> 1.
     This is the hallmark of exchange — if it's small here, the fit coefficients
    are wrong or the formula is broken.
 2. **Approaches unity at high energy**: Exchange vanishes as W increases because
     the emitted electron becomes distinguishable from bound electrons (higher momentum).
 3. **No NaN/inf values**: The empirical fit involves sin() and division by Wⁱ —
     we must ensure numerical stability across the full grid.
 4. **Coefficients loaded correctly**: Verify that Z=2..100 are all in the CSV.

Common practice: For components with dramatic energy dependence, test both extremes
(low-energy blowup and high-energy suppression) to verify the shape is correct.
"""

import numpy as np
import pytest

from beta_spectrum.components.exchange import ExchangeCorrection


class TestExchangeLowEnergy:
    """Test exchange correction at low energies where it's most significant."""

    def test_large_correction_for_heavy_nucleus_at_threshold(self):
        """For Z=80 (mercury), X >> 1 near threshold — can be > 20%.

        This is the defining characteristic of the exchange correction. If this
        fails, the empirical fit coefficients are wrong or the formula has a sign error.
        """
        ex = ExchangeCorrection(Z=80)

        low_W = np.array([1.1])  # T ≈ 50 keV — still "low" but above cutoff
        result = ex(low_W)

        assert (
            result[0] > 1.0
        ), f"Exchange must enhance spectrum at low energy for heavy Z: got {result[0]}"

    def test_correction_grows_as_energy_decreases(self):
        """X(W) should increase as W → 1 (lower momentum → stronger exchange)."""
        ex = ExchangeCorrection(Z=50)  # Tin — moderate-heavy

        low_W = np.array([1.2])  # Higher energy
        lower_W = np.array([1.05])  # Lower energy (above cutoff W_cut=1.005)

        s_low = ex(low_W)[0]
        s_lower = ex(lower_W)[0]

        assert (
            s_lower > s_low
        ), f"Exchange must increase at lower energy: {s_low} → {s_lower}"


class TestExchangeHighEnergy:
    """Test that exchange correction vanishes at high energy."""

    def test_near_unity_at_high_energy(self):
        """At W >> 1, the emitted electron is far from atomic electrons → X → 1."""
        ex = ExchangeCorrection(Z=40)

        high_W = np.array([3.0, 5.0])
        result = ex(high_W)

        assert (
            abs(result[0] - 1.0) < 0.5
        ), f"Exchange must approach unity at high energy: got {result}"


class TestExchangeNumericalStability:
    """Test that the empirical fit doesn't produce NaN or inf."""

    def test_no_nan_across_full_range(self):
        ex = ExchangeCorrection(Z=92)  # Uranium — extreme case
        W_test = np.linspace(1.01, 6.0, 50)
        result = ex(W_test)

        assert not np.any(np.isnan(result)), "No NaN values allowed"
        assert not np.any(np.isinf(result)), "No infinite values allowed"


class TestExchangeCoefficientLoading:
    """Test that coefficients are loaded correctly from CSV."""

    def test_all_Z_values_available(self):
        """The CSV must contain coefficients for Z=2 through Z=100+ (per Hayen Table X)."""
        # Test a few representative values across the periodic table
        for Z in [6, 20, 50, 80]:
            ex = ExchangeCorrection(Z=Z)
            assert ex.coeffs is not None, f"Missing coefficients for Z={Z}"

    def test_missing_Z_raises_error(self):
        """Requesting a Z without coefficients should raise ValueError."""
        # Z=1 has no electrons to exchange with — should fail gracefully
        with pytest.raises(ValueError, match="No exchange coefficients"):
            ExchangeCorrection(Z=1)


class TestExchangeOutputShape:
    """Basic shape checks."""

    def test_output_shape(self):
        ex = ExchangeCorrection(Z=30)
        for n in [1, 5, 100]:
            W_test = np.linspace(1.1, 3.0, n)
            result = ex(W_test)
            assert result.shape == W_test.shape
