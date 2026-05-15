"""
Shared API contract tests for all SpectrumComponent subclasses.

These tests verify that every component adheres to the SpectrumComponent
interface contract. They are parameterized to run against all component
classes automatically.

## Contract requirements

1. `__call__(W) -> np.ndarray` — the component is callable with W array
2. `__init__` accepts `logger: Optional[Logger] = None`
3. `__repr__` returns a non-empty string
4. Subclasses implement `__call__` (abstract method)

## Usage

Each component test file imports these helpers:

```python
import pytest
from beta_spectrum.components.fermi import FermiFunction
from tests.common.test_api_tests import test_component_is_callable

class TestFermiAPI:
    @pytest.mark.api
    def test_callable(self):
        ff = FermiFunction(Z=20, A=40)
        test_component_is_callable(ff, "FermiFunction")
```
"""

import numpy as np
import pytest


@pytest.mark.api
def _test_component_is_callable(component, component_name="Component"):
    """Component must be callable with array input."""
    assert callable(component), (
        f"{component_name}: component is not callable"
    )

    W = np.array([1.5, 2.0, 3.0])
    result = component(W)

    assert isinstance(result, np.ndarray), (
        f"{component_name}: __call__ returned {type(result)}, "
        f"expected np.ndarray"
    )
    assert result.shape == W.shape, (
        f"{component_name}: shape mismatch {W.shape} → {result.shape}"
    )


@pytest.mark.api
def _test_component_repr(component, component_name="Component"):
    """Component must have a non-empty repr."""
    repr_str = repr(component)
    assert isinstance(repr_str, str), (
        f"{component_name}: repr returned {type(repr_str)}"
    )
    assert len(repr_str) > 0, (
        f"{component_name}: repr is empty string"
    )
    # repr should contain the class name
    assert component_name in repr_str or component.__class__.__name__ in repr_str, (
        f"{component_name}: repr should contain class name. Got: {repr_str}"
    )


@pytest.mark.api
def _test_component_list_input(component, component_name="Component"):
    """Component must handle Python list input (numpy converts it)."""
    result = component([1.5, 2.0, 3.0])
    assert isinstance(result, np.ndarray), (
        f"{component_name}: list input returned {type(result)}"
    )


@pytest.mark.api
def _test_component_empty_input(component, component_name="Component"):
    """Component must handle empty array input gracefully."""
    result = component(np.array([]))
    assert isinstance(result, np.ndarray), (
        f"{component_name}: empty input returned {type(result)}"
    )
    assert result.shape == (0,), (
        f"{component_name}: empty input should return empty array"
    )
