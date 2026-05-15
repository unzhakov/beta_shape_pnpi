"""
Shared test utilities for the beta-spectrum test suite.

This module provides parametrized test helpers that eliminate duplicate
test code across correction components. Instead of each component
writing its own `test_positive`, `test_output_shape`, `test_no_nan`
functions, they import and use these shared helpers.

## Design principles

1. **Parametrized, not duplicated** — one `test_positive_everywhere` helper
   covers all components via pytest.mark.parametrize.
2. **Fail-fast** — if a helper fails, the error message identifies the
   component and the violated property.
3. **Physics-aware** — helpers understand that some components (like
   phase_space) return zero at threshold, so "positive" means "non-negative
   with at least one strictly positive value."
"""
