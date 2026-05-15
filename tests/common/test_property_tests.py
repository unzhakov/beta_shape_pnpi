"""
Shared property-based tests for all SpectrumComponent subclasses.

These tests verify universal properties that every correction component
must satisfy. Instead of duplicating `test_positive`, `test_output_shape`,
`test_no_nan` in every component test file, each file imports and
parametrizes these helpers.

## Property test categories

1. **Positivity** — multiplicative factors must be ≥ 0 (strictly > 0 except
   at physical boundaries like W=1 or W=W₀)
2. **Shape preservation** — output shape must match input shape
3. **Numerical stability** — no NaN or inf values across the full energy range
4. **Type correctness** — output is always numpy.ndarray

## Usage example

```python
import pytest
from beta_spectrum.components.fermi import FermiFunction
from tests.common.test_property_tests import (
    test_positive_values,
    test_output_shape_preserved,
    test_no_nan_across_range,
)

class TestFermiProperties:
    @pytest.mark.parametrize("W", [[1.5], [1.1, 2.0, 3.0], [1.01, 2.0, 3.0, 4.0, 5.0]])
    @pytest.mark.property
    def test_positive(self, W):
        ff = FermiFunction(Z=20, A=40)
        test_positive_values(ff, W, component_name="FermiFunction")

    @pytest.mark.parametrize("n", [1, 5, 100])
    @pytest.mark.property
    def test_shape(self, n):
        ff = FermiFunction(Z=20, A=40)
        W = np.linspace(1.1, 3.0, n)
        test_output_shape_preserved(ff, W, component_name="FermiFunction")
```

## Test markers

- `@pytest.mark.property` — marks a test as checking a universal property
- `@pytest.mark.quality` — marks a test as checking code quality (shape, type)

These markers allow running only property tests (`pytest -m property`)
or only quality tests (`pytest -m quality`) for fast feedback during development.
"""

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Property tests: positivity
# ---------------------------------------------------------------------------

@pytest.mark.property
def _test_positive_values(component, W, component_name="Component"):
    """Component must return positive values at all W in the physical range.

    For multiplicative correction factors, negative values indicate a
    fundamental bug (wrong sign in formula, overflow to negative, etc.).

    Parameters
    ----------
    component : SpectrumComponent or callable
        The component to test. Must be callable as `component(W)`.
    W : array_like
        Energy grid to evaluate at. Values should be in the physical range
        (W > 1 for β⁻, W < W₀ for endpoint).
    component_name : str
        Display name for error messages.
    """
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    assert not np.any(np.isnan(result)), (
        f"{component_name}: NaN values found in output. "
        f"W range: [{W_arr.min():.3f}, {W_arr.max():.3f}], "
        f"result range: [{np.nanmin(result):.6f}, {np.nanmax(result):.6f}]"
    )
    assert not np.any(np.isinf(result)), (
        f"{component_name}: Inf values found in output."
    )
    assert np.all(result > 0), (
        f"{component_name}: non-positive values found. "
        f"Min = {np.min(result):.10f} at W = {W_arr[np.argmin(result)][0]:.4f}"
    )


@pytest.mark.property
def _test_non_negative_values(component, W, component_name="Component", allow_zero=True):
    """Component must return non-negative values (allows zero at physical boundaries).

    Use this for components that can legitimately be zero (e.g., phase_space
    at W=1 or W=W₀).

    Parameters
    ----------
    allow_zero : bool
        If True, zero values are allowed. If False, strictly positive is required.
    """
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    assert not np.any(np.isnan(result)), (
        f"{component_name}: NaN values found in output."
    )
    assert not np.any(np.isinf(result)), (
        f"{component_name}: Inf values found in output."
    )

    if allow_zero:
        assert np.all(result >= 0), (
            f"{component_name}: negative values found. "
            f"Min = {np.min(result):.10f}"
        )
    else:
        assert np.all(result > 0), (
            f"{component_name}: non-positive values found. "
            f"Min = {np.min(result):.10f}"
        )


# ---------------------------------------------------------------------------
# Property tests: shape preservation
# ---------------------------------------------------------------------------

@pytest.mark.quality
def _test_output_shape_preserved(component, W, component_name="Component"):
    """Output shape must match input shape — no flattening or reshaping."""
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    assert isinstance(result, np.ndarray), (
        f"{component_name}: output is {type(result)}, expected np.ndarray"
    )
    assert result.shape == W_arr.shape, (
        f"{component_name}: shape mismatch. "
        f"Input {W_arr.shape} → Output {result.shape}"
    )


@pytest.mark.quality
def _test_scalar_input_returns_array(component, W_scalar, component_name="Component"):
    """Single-element array input should return array (not scalar)."""
    result = component(W_scalar)
    assert isinstance(result, np.ndarray), (
        f"{component_name}: single-element input returned {type(result)}, "
        f"expected np.ndarray"
    )
    assert result.ndim == 1, (
        f"{component_name}: expected 1D output, got {result.ndim}D"
    )


# ---------------------------------------------------------------------------
# Property tests: numerical stability
# ---------------------------------------------------------------------------

@pytest.mark.quality
def _test_no_nan_across_range(component, W, component_name="Component"):
    """No NaN values across the entire energy range."""
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    nan_mask = np.isnan(result)
    if np.any(nan_mask):
        nan_W = W_arr[nan_mask]
        pytest.fail(
            f"{component_name}: NaN values at W = {nan_W}. "
            f"Result min = {np.nanmin(result):.6f}, max = {np.nanmax(result):.6f}"
        )


@pytest.mark.quality
def _test_no_inf_across_range(component, W, component_name="Component"):
    """No infinite values across the entire energy range."""
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    inf_mask = np.isinf(result)
    if np.any(inf_mask):
        inf_W = W_arr[inf_mask]
        pytest.fail(
            f"{component_name}: Inf values at W = {inf_W}. "
            f"Result min = {np.nanmin(result):.6f}, max = {np.nanmax(result):.6f}"
        )


@pytest.mark.quality
def _test_full_range_stability(component, W, component_name="Component"):
    """Full energy range: no NaN, no inf, finite values."""
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    assert np.all(np.isfinite(result)), (
        f"{component_name}: non-finite values in full range. "
        f"NaN count: {np.sum(np.isnan(result))}, "
        f"Inf count: {np.sum(np.isinf(result))}"
    )


# ---------------------------------------------------------------------------
# Property tests: type correctness
# ---------------------------------------------------------------------------

@pytest.mark.quality
def _test_output_is_ndarray(component, W, component_name="Component"):
    """Output must be numpy.ndarray for any valid input."""
    W_arr = np.asarray(W, dtype=float)
    result = component(W_arr)

    assert isinstance(result, np.ndarray), (
        f"{component_name}: output type is {type(result).__name__}, "
        f"expected np.ndarray"
    )


# ---------------------------------------------------------------------------
# Parametrized test templates for quick component integration
# ---------------------------------------------------------------------------

def make_positive_test(component_class, test_params, component_name=None):
    """Factory function to create a parametrized positivity test.

    Parameters
    ----------
    component_class : type
        The component class to instantiate (e.g., FermiFunction).
    test_params : list[dict]
        List of keyword arguments for component construction.
        Each dict is tested at multiple W values.
    component_name : str, optional
        Display name. Defaults to component_class.__name__.

    Example
    -------
    tests = make_positive_test(
        FermiFunction,
        [{"Z": 1, "A": 1}, {"Z": 20, "A": 40}, {"Z": 92, "A": 238}],
    )
    """
    if component_name is None:
        component_name = component_class.__name__

    W_values = [
        np.array([1.5]),
        np.linspace(1.1, 3.0, 5),
        np.linspace(1.01, 5.0, 20),
    ]

    @pytest.mark.property
    @pytest.mark.parametrize("params", test_params)
    @pytest.mark.parametrize("W", W_values)
    def test_positive(component_params, W_arr):
        comp = component_class(**component_params)
        test_positive_values(comp, W_arr, component_name)

    return test_positive


def make_stability_test(component_class, test_params, component_name=None):
    """Factory function to create a parametrized stability test."""
    if component_name is None:
        component_name = component_class.__name__

    W_values = [
        np.linspace(1.01, 6.0, 50),
        np.linspace(1.001, 10.0, 100),
    ]

    @pytest.mark.quality
    @pytest.mark.parametrize("params", test_params)
    @pytest.mark.parametrize("W", W_values)
    def test_stability(component_params, W_arr):
        comp = component_class(**component_params)
        test_full_range_stability(comp, W_arr, component_name)

    return test_stability
