---
title: "Test Architecture — Organization & Strategy"
date: 2026-05-15
tags:
  - testing/architecture
status: active
aliases: [test strategy, test categories, testing]
---

# Test Architecture — How Tests Are Organized

## Overview

The test suite currently has **230 tests** across 14 files. The tests fall into two fundamentally different categories:

### Category 1: Code Quality Tests (CQT)

Tests that verify the code works as a correct Python program:
- **Signatures**: function accepts correct arguments, correct types
- **Return types**: output is `np.ndarray`, correct shape
- **Error handling**: `ValueError` raised for invalid inputs
- **API contracts**: factory methods work, config parsing works

These tests are **implementation-agnostic** — they don't care about physics, just that the code doesn't crash and returns the right types.

### Category 2: Physics Tests (PT)

Tests that verify the code implements correct physics:
- **Physical properties**: positivity, monotonicity, limits
- **Numerical stability**: no NaN/inf at boundaries
- **Known limits**: Z=1 → F₀=1, light nucleus → correction ≈ 1
- **Energy dependence**: screening increases at low energy, radiative increases toward endpoint
- **Integration**: full pipeline produces valid spectrum

## Current State Analysis

### Redundant Patterns (18 instances)

Every correction component has:
1. `test_positive*` — checks output > 0
2. `test_output_shape*` — checks output shape matches input
3. `test_no_nan*` or `test_stability*` — checks no NaN/inf

These are **mechanically identical** across `test_fermi.py`, `test_finite_size.py`, `test_screening.py`, `test_exchange.py`, `test_radiative.py`, `test_phase_space.py`.

### Test Organization (Current)

```
tests/
├── conftest.py                 ← shared fixtures (W_low, W_mid, W_high, W_full, large_W)
├── test_fermi.py               ← 7 tests — Fermi function
├── test_finite_size.py         ← 9 tests — L0, U corrections
├── test_screening.py           ← 5 tests — atomic screening
├── test_exchange.py            ← 7 tests — atomic exchange
├── test_radiative.py           ← 17 tests — radiative correction (largest)
├── test_phase_space.py         ← 8 tests — baseline spectrum
├── test_detector_response.py   ← 23 tests — detector models
├── test_spectrum.py            ← 14 tests — full pipeline integration
├── test_fitter.py              ← 16 tests — χ² minimization
├── test_cw_extractor.py        ← 18 tests — shape factor extraction
├── test_multi_branch.py        ← 26 tests — multi-branch decay support
├── test_nuclear_data.py        ← 21 tests — ENSDF/JSON input
├── test_logging_utils.py       ← 23 tests — logging infrastructure
└── (empty)
```

### Test Distribution by Type

| File | CQT (code quality) | PT (physics) | Total |
|---|---|---|---|
| test_fermi.py | 3 (shape, error) | 4 (positive, stability, monotonicity) | 7 |
| test_finite_size.py | 2 (shape) | 7 (positive, limits, Z-scaling) | 9 |
| test_screening.py | 1 (shape) | 4 (positive, energy dep, limits) | 5 |
| test_exchange.py | 2 (shape, error) | 5 (positive, stability, energy dep) | 7 |
| test_radiative.py | 2 (shape) | 15 (positive, stability, endpoint, Z-dep) | 17 |
| test_phase_space.py | 2 (shape, type) | 6 (positive, limits, neutrino mass) | 8 |
| test_detector_response.py | 5 (shape, error, repr) | 18 (normalization, peak, resolution) | 23 |
| test_spectrum.py | 6 (API, config, shape) | 8 (pipeline, integration) | 14 |
| test_fitter.py | 4 (API, config) | 12 (fit quality, covariance, bounds) | 16 |
| test_cw_extractor.py | 4 (API, config) | 14 (extraction accuracy, Kurie, gV/gA) | 18 |
| test_multi_branch.py | 3 (API, config) | 23 (physics, branching, normalization) | 26 |
| test_nuclear_data.py | 6 (API, error handling) | 15 (parsing, conversion, integration) | 21 |
| test_logging_utils.py | 23 | 0 | 23 |
| **Total** | **79** | **151** | **230** |

## Proposed New Structure

### Directory Layout

```
tests/
├── conftest.py                 ← shared fixtures
├── pytest.ini                  ← custom markers (see below)
├── common/                     ← shared test utilities
│   ├── __init__.py
│   ├── test_property_tests.py  ← property-based test helpers (positive, shape, NaN)
│   └── test_api_tests.py       ← API contract test helpers
├── physics/                    ← Physics tests (PT)
│   ├── __init__.py
│   ├── test_phase_space.py     ← Phase space physics
│   ├── test_fermi.py           ← Fermi function physics
│   ├── test_finite_size.py     ← L0, U physics
│   ├── test_screening.py       ← Screening physics
│   ├── test_exchange.py        ← Exchange physics
│   ├── test_radiative.py       ← Radiative physics
│   ├── test_detector_response.py ← Detector physics (normalization, resolution)
│   ├── test_spectrum.py        ← Full pipeline physics
│   ├── test_fitter.py          ← Fitter physics (fit quality)
│   ├── test_cw_extractor.py    ← Shape factor extraction physics
│   └── test_multi_branch.py    ← Multi-branch physics
├── quality/                    ← Code Quality tests (CQT)
│   ├── __init__.py
│   ├── test_fermi.py           ← Fermi API/signatures
│   ├── test_finite_size.py     ← L0, U API
│   ├── test_screening.py       ← Screening API
│   ├── test_exchange.py        ← Exchange API
│   ├── test_radiative.py       ← Radiative API
│   ├── test_phase_space.py     ← Phase space API
│   ├── test_detector_response.py ← Detector API
│   ├── test_spectrum.py        ← BetaSpectrum API
│   ├── test_fitter.py          ← CurveFitter API
│   ├── test_cw_extractor.py    ← CWExtractor API
│   ├── test_multi_branch.py    ← Multi-branch API
│   ├── test_nuclear_data.py    ← ENSDF/JSON API
│   └── test_logging_utils.py   ← Logging API
└── integration/                ← Cross-component tests
    ├── __init__.py
    └── test_full_pipeline.py   ← End-to-end from config to CSV export
```

### Pytest Markers (in pyproject.toml)

```toml
[tool.pytest.ini_options]
markers = [
    "property: Physics property tests (positivity, monotonicity, limits)",
    "physics: Physics behavior tests (energy dependence, Z dependence)",
    "api: API contract tests (signatures, return types, error handling)",
    "quality: Code quality tests (shape, type, NaN, stability)",
    "integration: Full pipeline integration tests",
    "slow: Tests that take >1s to run",
    "requires_ensdf: Tests requiring paceENSDF package",
]
```

### Key Design Decisions

1. **Property tests go in `common/test_property_tests.py`** — a shared module with parametrized test functions. Each correction component imports and uses these instead of writing duplicate `test_positive`, `test_output_shape`, `test_no_nan` tests.

2. **Physics tests stay in `physics/`** — each file tests the unique physical behavior of that component.

3. **API tests go in `quality/`** — tests that check function signatures, return types, error handling. These are boring but essential.

4. **Integration tests in `integration/`** — tests that span multiple components (e.g., full pipeline from config to CSV export).

## Implementation Plan

### Phase 1: Create infrastructure (no behavior change)
1. Create `tests/common/`, `tests/physics/`, `tests/quality/`, `tests/integration/`
2. Add `pytest.ini` with custom markers
3. Move `conftest.py` to shared location
4. Create `common/test_property_tests.py` with parametrized property tests

### Phase 2: Migrate physics tests
1. Move all physics tests from current files to `tests/physics/`
2. Replace duplicate `test_positive`, `test_output_shape`, `test_no_nan` with calls to `common/test_property_tests.py` helpers
3. Keep unique physics tests in place

### Phase 3: Migrate API/quality tests
1. Move all API/quality tests to `tests/quality/`
2. Replace with parametrized API contract tests

### Phase 4: Create integration tests
1. Create `tests/integration/test_full_pipeline.py`
2. Test: config → spectrum → convolve → export → reimport
3. Test: multi-branch → weighted sum → export

### Phase 5: Cleanup
1. Remove old test files
2. Update `pyproject.toml` pytest config
3. Verify all 230 tests still pass
4. Document the new structure in `docs/test-architecture.md`

## References

- [pytest documentation on markers](https://docs.pytest.org/en/stable/how-to/mark.html)
- [NumPy testing guidelines](https://numpy.org/doc/stable/reference/testing.html)
- [Property-based testing with Hypothesis](https://hypothesis.readthedocs.io/)
- [SciPy test organization](https://docs.scipy.org/doc/scipy/reference/testing.html)
