---
title: "Test Migration TODO — From Flat to Structured Suite"
date: 2026-05-15
tags:
  - testing/todo
  - testing/migration
status: active
aliases: [test migration, test reorganization]
related_notes:
  - test-architecture
---

# Test Migration TODO

> [!info] Current state
> - **Old tests**: 230 tests in `tests/*.py` (flat structure)
> - **New tests**: 94 tests in `tests/{physics,quality,common}/` (structured)
> - **New tests already passing**: ✅ 94/94
> - **Old tests still passing**: need to verify after removing duplicates

## Phase 1: Remove duplicate old tests (P0 — blocking)

The old files in `tests/` root contain tests that are now covered by the new structure. Remove them:

- [x] **1.1** Remove `tests/test_fermi.py` — covered by `tests/physics/test_fermi.py` + `tests/quality/test_api_fermi.py`
- [x] **1.2** Remove `tests/test_finite_size.py` — covered by `tests/physics/test_finite_size.py` + `tests/quality/test_api_finite_size.py`
- [x] **1.3** Remove `tests/test_screening.py` — covered by `tests/physics/test_screening.py` + `tests/quality/test_api_screening.py`
- [x] **1.4** Remove `tests/test_exchange.py` — covered by `tests/physics/test_exchange.py` + `tests/quality/test_api_exchange.py`
- [x] **1.5** Remove `tests/test_radiative.py` — covered by `tests/physics/test_radiative.py` + `tests/quality/test_api_radiative.py`
- [x] **1.6** Remove `tests/test_phase_space.py` — covered by `tests/physics/test_phase_space.py` + `tests/quality/test_api_phase_space.py`

**After removal**: ~136 tests remain in new structure (fitter, spectrum, detector, cw_extractor, multi_branch, nuclear_data, logging).

## Phase 2: Migrate remaining old tests (P1 — high priority)

Move the remaining old tests to the new structure:

### 2.1 Integration tests → `tests/integration/`
- [x] **2.1.1** `tests/test_spectrum.py` → `tests/integration/test_spectrum_pipeline.py`
  - Full pipeline tests (from_config → evaluate → components)
  - Declarative detector response tests
  - These are integration tests, not component tests

### 2.2 Physics tests → `tests/physics/`
- [x] **2.2.1** `tests/test_detector_response.py` → `tests/physics/test_detector_response.py` (split into physics + quality)
  - 23 tests — physics behavior (normalization, resolution, convolution)
  - Split into: `test_detector_response.py` (physics) + `tests/quality/test_api_detector_response.py` (API)

- [x] **2.2.2** `tests/test_fitter.py` → `tests/physics/test_fitter.py`
  - 16 tests — fit quality, covariance, confidence intervals
  - Split: physics (fit quality) + quality (API)

- [x] **2.2.3** `tests/test_cw_extractor.py` → `tests/physics/test_cw_extractor.py`
  - 18 tests — shape factor extraction accuracy
  - Split: physics (extraction) + quality (API)

- [x] **2.2.4** `tests/test_multi_branch.py` → `tests/physics/test_multi_branch.py`
  - 26 tests — multi-branch physics
  - Split: physics (branching) + quality (API)

- [x] **2.2.5** `tests/test_nuclear_data.py` → `tests/quality/test_nuclear_data.py`
  - 21 tests — mostly API/parsing tests
  - These are quality tests (input validation, parsing)

### 2.3 Quality tests → `tests/quality/`
- [x] **2.3.1** `tests/test_logging_utils.py` → `tests/quality/test_logging_utils.py`
  - 23 tests — all API tests (no physics)
  - Direct move, no splitting needed

## Phase 3: Create shared integration test (P2)
- [ ] **3.1** Create `tests/integration/test_full_pipeline.py`
  - End-to-end: config → spectrum → convolve → export → reimport
  - Multi-branch → weighted sum → export
  - Verify CSV export matches reimport

## Phase 4: Clean up (P3)
- [x] **4.1** Remove old `tests/test_*.py` files that were migrated
- [x] **4.2** Verify all tests pass: `python3 -m pytest tests/ -v` ✅ 256 passed
- [x] **4.3** Verify test count matches expected
- [x] **4.4** Update `docs/test-architecture.md` with final structure
- [x] **4.5** Add `.gitkeep` to `tests/common/`, `tests/integration/` directories

## Phase 5: Advanced improvements (P4 — nice to have)
- [x] **5.1** Add `pytest-timeout` (30s, thread-based) to pyproject.toml
- [x] **5.2** Add `pytest.mark.requires_ensdf` marker to ENSDF tests
- [x] **5.3** Add `hypothesis` for property-based edge case discovery
  - 5 new tests in `tests/common/test_hypothesis_tests.py`
  - Random Z, A, W for Fermi function
  - Random W0, W for phase space
  - Hypothesis already found real bugs (see commit message)
- [ ] **5.4** Add `pytest-timeout` to prevent hangs (already done in 5.1)
- [ ] **5.5** Add CI badge for test status (nice to have)
- [ ] **5.1** Add `pytest.mark.slow` to long-running tests
- [ ] **5.2** Add `pytest.mark.requires_ensdf` to tests needing paceENSDF
- [ ] **5.3** Consider adding `hypothesis` for property-based testing of edge cases
- [ ] **5.4** Add `pytest-timeout` to prevent hangs
- [ ] **5.5** Add CI badge for test status

## Files to create/modify

### New files (already created)
- ✅ `tests/common/__init__.py`
- ✅ `tests/common/test_property_tests.py` — shared property helpers
- ✅ `tests/common/test_api_tests.py` — shared API helpers
- ✅ `tests/physics/test_fermi.py`
- ✅ `tests/physics/test_phase_space.py`
- ✅ `tests/physics/test_finite_size.py`
- ✅ `tests/physics/test_screening.py`
- ✅ `tests/physics/test_exchange.py`
- ✅ `tests/physics/test_radiative.py`
- ✅ `tests/quality/test_api_fermi.py`
- ✅ `tests/quality/test_api_phase_space.py`
- ✅ `tests/quality/test_api_finite_size.py`
- ✅ `tests/quality/test_api_screening.py`
- ✅ `tests/quality/test_api_exchange.py`
- ✅ `tests/quality/test_api_radiative.py`

### Files to create (remaining)
- `tests/integration/test_spectrum_pipeline.py`
- `tests/integration/test_full_pipeline.py`
- `tests/physics/test_detector_response.py`
- `tests/quality/test_api_detector_response.py`
- `tests/physics/test_fitter.py`
- `tests/quality/test_api_fitter.py`
- `tests/physics/test_cw_extractor.py`
- `tests/quality/test_api_cw_extractor.py`
- `tests/physics/test_multi_branch.py`
- `tests/quality/test_api_multi_branch.py`
- `tests/quality/test_nuclear_data.py`
- `tests/quality/test_logging_utils.py`

### Files to delete
- `tests/test_fermi.py`
- `tests/test_finite_size.py`
- `tests/test_screening.py`
- `tests/test_exchange.py`
- `tests/test_radiative.py`
- `tests/test_phase_space.py`
- `tests/test_spectrum.py`
- `tests/test_detector_response.py`
- `tests/test_fitter.py`
- `tests/test_cw_extractor.py`
- `tests/test_multi_branch.py`
- `tests/test_nuclear_data.py`
- `tests/test_logging_utils.py`

## Test count comparison

| Category | Before | After (target) |
|---|---|---|
| Common helpers | 0 | ~12 (shared, not counted) |
| Physics tests | ~151 | ~151 (same, reorganized) |
| Quality/API tests | ~79 | ~79 (same, reorganized) |
| Integration tests | ~0 | ~10 (new) |
| **Total unique tests** | **230** | **~240** |

The duplicate tests (positivity, shape, NaN across 6 components = ~18 tests) are now covered by shared helpers, so the **effective test count decreases** while coverage stays the same or improves.

## References

- [pytest markers documentation](https://docs.pytest.org/en/stable/how-to/mark.html)
- [pytest test collection](https://docs.pytest.org/en/stable/explanation/goodpractices.html#test-collection)
- [NumPy testing guidelines](https://numpy.org/doc/stable/reference/testing.html)
- [Property-based testing with Hypothesis](https://hypothesis.readthedocs.io/)
