# Strategic Development Roadmap

## Scientific Goal

Extract the shape factor C(W) from the experimentally measured beta spectrum of ⁹⁹Tc, and parametrize C(W) via the vector coupling constant g_V from the theory of beta decay.

This requires repeating the analysis pipeline of [Paulsen et al., Phys. Rev. C **110**, 055503 (2024)](https://doi.org/10.1103/PhysRevC/110/055503), which used Metallic Magnetic Calorimeters (MMCs) to measure the ⁹⁹Tc beta spectrum with sub-keV threshold and ~100 eV resolution.

See `docs/refs/2024_Paulsen_T99_beta_spectrum/summary.md` for full analysis details.

______________________________________________________________________

## Package Architecture (v0.4.0 — In Progress)

### Current State (v0.3.x)

Everything lives in `beta_spectrum/` — a monolithic package where theory, fitting, detector response, and analysis are intermingled. This is becoming unwieldy.

### Target Architecture (v0.4.0)

```
beta_spectrum/              ← theory-only core
├── __init__.py
├── base.py                 ← SpectrumComponent base
├── constants.py            ← physical constants
├── utils.py                ← T_to_W, W_to_T, momentum
├── spectrum.py             ← BetaSpectrum, SpectrumConfig, BranchConfig
├── components/             ← individual correction factors
│   ├── phase_space.py
│   ├── fermi.py
│   ├── finite_size.py
│   ├── screening.py
│   ├── exchange.py
│   ├── radiative.py
│   └── detector_response.py   ← stays here (part of theory)
└── nuclear_data.py         ← paceENSDF integration

exp_data/                   ← experimental data handling
├── __init__.py
├── raw_data.py             ← raw data loading (ROOT, CSV, ASCII)
├── calibration.py          ← energy calibration, line fitting
├── fitters.py              ← simple Gaussian fitters for cal sources
├── corrections.py          ← dead-time, pile-up, background
└── spectrum.py             ← experimental spectrum container

fitter/                     ← main analysis framework
├── __init__.py
├── model.py                ← theoretical model + detector convolution
├── fit_engine.py           ← multi-parameter fitting orchestration
├── extractor.py            ← C(W), g_V, g_A extraction
├── result.py               ← FitResult, analysis results
└── report.py               ← PDF report generation

analyzer/                   ← visualization (stays as-is)
└── analyzer.py             ← BetaSpectrumAnalyzer
```

### Rationale

- **`beta_spectrum`** — pure theory. Given a nuclide and corrections, produces a theoretical spectrum. No fitting, no data, no detector convolution.
- **`exp_data`** — experimental data handling. Loading raw detector outputs, calibrating energy scale, fitting calibration peaks, applying instrumental corrections.
- **`fitter`** — the analysis workhorse. Takes theory model + detector response + experimental data → performs convolution + fitting → extracts physics parameters.

### Migration Plan

#### Phase 1: Create `exp_data` module
- [ ] `exp_data/raw_data.py` — load experimental spectra from various formats (CSV, ROOT, ASCII)
- [ ] `exp_data/calibration.py` — energy calibration using known lines (X-rays, gammas)
- [ ] `exp_data/fitters.py` — simple Gaussian + background fitters for calibration peaks
- [ ] `exp_data/corrections.py` — dead-time, pile-up, background subtraction
- [ ] `exp_data/spectrum.py` — `ExpSpectrum` dataclass: energies, counts, errors, metadata
- [ ] Move experimental data tests from `tests/integration/` → `tests/exp_data/`

#### Phase 2: Create `fitter` module
- [ ] `fitter/model.py` — `TheoreticalModel` class: wraps `BetaSpectrum` + `DetectorResponse`, provides `evaluate(E)` that includes convolution
- [ ] `fitter/fit_engine.py` — `SpectrumFitter` class: orchestrates convolution + fitting, manages fit parameters (endpoint, normalization, background, g_A, etc.)
- [ ] `fitter/extractor.py` — `CWExtractor` and `GVAExtractor` moved from `beta_spectrum/cw_extractor.py`
- [ ] `fitter/result.py` — `AnalysisResult` class: wraps `FitResult` with physics-specific metadata (nuclide, endpoint, extracted C(W), g_V, g_A)
- [ ] Move `CurveFitter` / `FitConfig` / `FitResult` from `beta_spectrum/fitter.py` → `fitter/fit_engine.py` (or keep as shared utility)
- [ ] Move detector response convolution from `beta_spectrum/components/detector_response.py` → `fitter/model.py` (or keep `DetectorResponse` class in `beta_spectrum` as a data structure)

#### Phase 3: Clean up `beta_spectrum`
- [ ] Remove `cw_extractor.py` (moved to `fitter/`)
- [ ] Remove `fitter.py` (moved to `fitter/`)
- [ ] Move `DetectorResponse` to `beta_spectrum/components/detector_response.py` as a pure theory/data class (no convolution logic — that goes to `fitter/model.py`)
- [ ] Clean up `__init__.py` — only export theory classes + DetectorResponse data class
- [ ] Remove CLI `bs_pnpi` from `beta_spectrum` (move to a `bin/` or `scripts/` directory as a standalone tool)

#### Phase 4: Update package structure
- [ ] Update `pyproject.toml` — add `exp_data` and `fitter` as subpackages
- [ ] Update `__init__.py` at root level if needed
- [ ] Move CLI entry point to a separate package or keep as `beta_spectrum.cli` with re-exports
- [ ] Update all imports across tests
- [ ] Run full test suite
- [ ] Update documentation

______________________________________________________________________

## Track A: ⁹⁹Tc Spectrum-Shape Analysis

### A1. Research & Data Collection

- [ ] Parse additional PDF sources on ⁹⁹Tc beta decay and spectrum-shape method
  - Kostensalo & Suhonen, Phys. Rev. C **96**, 024317 (2017) — theoretical spectrum shape predictions
  - Suhonen, Frontier. Phys. **5**, 55 (2017) — review of g_A quenching
  - Haaranen et al., Phys. Rev. C **93**, 034308 (2016) — effective coupling constants
  - Additional experimental ⁹⁹Tc spectra (historical data for comparison)
- [ ] Research methods to extract g_V from fitted C(W) shape
  - Study the relationship between C(W) parametrization and g_V, g_A effective values
  - Understand how CVC (Conserved Vector Current) constrains g_V
  - Review the spectrum-shape method formalism (Behrens-Bühring + shell model)
- [ ] Document the full analysis pipeline: experimental data → theoretical model → C(W) extraction → g_V extraction

### A2. Detector Response Function

- [x] ~~Implement analytical detector response function for 4π semiconductor detector~~
  - ~~Gaussian core with low-energy tail (charge collection effects)~~
  - ~~Parameters: energy resolution σ(E), tail fraction, tail shape~~
  - ~~Energy-dependent resolution σ(E) = a + b·√E~~
  - ~~Implemented: Gaussian, Gaussian+tail, Tikhonov models with energy-dependent resolution and Fano factor~~
- [x] ~~Implement convolution/deconvolution routines~~
  - ~~Convolve theoretical spectrum with detector response for comparison to data~~
  - ~~Deconvolve measured spectrum to recover true shape (for initial C(W) extraction)~~
  - ~~Iterative unfolding (Richardson-Lucy or similar)~~
  - ~~Implemented: convolve() and convolve_batch() with tabulated response support~~
- [ ] Monte Carlo simulation (GEANT4) — out of scope for this project but noted for future
  - Full detector geometry simulation
  - Energy loss in source, dead layers, absorber
  - Response matrix for unfolding (as done in Paulsen et al. with EGSNRC)

### A3. Data Processing & Analysis

- [ ] Data processing pipeline for 1-hour run spectra
  - Background subtraction
  - Energy calibration (using known X-ray/gamma lines)
  - Dead-time correction
  - Pulse pile-up correction
- [ ] Endpoint fitting and energy calibration
  - Kurie plot analysis for endpoint determination
  - Polynomial calibration using known spectral lines
  - Linearity correction across energy range
- [ ] Data quality assessment procedures
  - Residual analysis (data vs model)
  - χ²/ndf evaluation
  - Consistency checks across multiple runs
  - Statistical combination of long-exposure spectra

### A4. Fitter for C(W) Extraction

- [x] ~~Implement fitter routine to extract experimental C(W) from data~~
  - ~~χ² minimization: theoretical spectrum (with all corrections) vs measured spectrum~~
  - ~~Free parameters: C(W) shape parameters, endpoint energy, normalization, background~~
  - ~~Covariance matrix and uncertainty propagation~~
  - ~~Implemented: `CurveFitter` with `least_squares` optimization, `FitResult` with covariance, chi2, residuals~~
- [x] ~~Parametrize C(W) in terms of g_V (and g_A)~~
  - ~~Fit C(W) data to theoretical parametrization~~
  - ~~Extract g_V^eff and g_A^eff from the fit~~
  - ~~Compare with Paulsen et al. results: g_A^eff = 0.574(36), g_V^eff = 0.376(5)~~
  - ~~Implemented: `CWExtractor` with Kurie plot analysis, parametrized fitting, g_V/g_A extraction~~
- [ ] Systematic uncertainty analysis
  - Vary correction implementations within uncertainties
  - Test sensitivity to detector response model
  - Background model dependence

______________________________________________________________________

## Track B: Package Usability Improvements

### B1. Nuclear Data Integration

- [x] ~~paceENSDF integration for automated nuclear data retrieval~~
  - ~~`get_decay_info_from_paceENSDF(nuclide, decay_type)` — retrieves Q-value, half-life, spin/parity, branches, forbiddenness from ENSDF~~
  - ~~`decay_info_to_config()` — converts DecayInfo to SpectrumConfig with branch support~~
  - ~~`create_config_from_source("paceENSDF", nuclide="Tc99")` — one-liner config creation~~
  - ~~`_parse_nuclide_symbol()` — parses "Tc99" → (Tc, 43, 99) with full Z lookup for Z=1..98~~
  - ~~`_resolve_decay_index()` — finds correct decay index (ground state by default)~~
  - ~~`_get_decay_mode_symbol()` — maps beta_minus→BM, beta_plus/ec→ECBP~~
  - ~~FORBIDDENNESS_MAP — translates paceENSDF codes (0A, 1F, 2F, …) to transition_type~~
  - ~~CLI `--nuclide Tc99` uses paceENSDF automatically~~
  - ~~Branch intensities, log_ft, level energies, transition types all populated from ENSDF~~
  - ~~Tested with Tc99 (2 branches, F2) and Co60 (multi-branch)~~

### B2. Input Flexibility

- [x] ~~Generalize parameter inputs via custom input file (YAML/JSON)~~
  - ~~Declare isotope, transition type, detector parameters~~
  - ~~Toggle corrections on/off per-component~~
  - ~~Reproducible calculation configurations~~
  - ~~Implemented: `load_json_input()`, `json_to_config()`, `DEFAULT_JSON_SCHEMA` with full detector param support, CLI `--input` flag~~
- [x] ~~Retrieve parameters directly from ENSDF database~~
  - ~~Auto-populate W0, Z, A, transition type from isotope name~~
  - ~~Reduce manual configuration errors~~
  - ~~Implemented: `get_decay_info_from_paceENSDF()` returns DecayInfo with Z_parent, Z_daughter, A_number, endpoint_MeV, transition_type, forbiddenness_code, half_life, branches; CLI `--nuclide` uses this automatically~~

### B3. Multiple Branch Decays

- [x] ~~Add support for multiple branch decays~~
  - ~~Calculate total spectral shape weighted by branch intensities~~
  - ~~Handle decays with significant branching to excited states~~
  - ~~Sum contributions from all branches with proper weighting~~
  - ~~Implemented: `BranchConfig` dataclass, multi-branch `BetaSpectrum` with per-branch calculators~~
  - ~~Branch data from paceENSDF: intensity, endpoint (Q-value - excitation energy), transition type~~
  - ~~Intensity normalization to 100% across all branches~~
  - ~~CSV export: total spectrum, per-branch spectra, per-branch components (PhaseSpace, Radiative), universal components (Fermi, Screening, etc.)~~
  - ~~Plot: total spectrum + individual branches with legend (transition type, endpoint energy)~~
  - ~~Debug plot: vertical layout with per-branch panels, universal components panel~~
  - ~~Energy grid extends to max branch endpoint, per-branch spectra masked beyond their endpoint~~
  - ~~Single-branch = trivial case of multi-branch (no separate mode)~~
  - ~~Per-branch transition types extracted from ENSDF (not defaulting to 'A')~~
  - ~~`--intensity-cutoff` CLI option to filter out negligible branches~~
  - ~~Branch spectra plotted as raw (un-normalized) contributions so total = sum of branches~~

### B4. Beta+ Decay Support

- [ ] Modify corrections to support β⁺ decays
  - Fermi function sign change for positrons
  - Radiative corrections for β⁺ (Sirlin corrections differ for e⁻ vs e⁺)
  - Screening correction for positrons (different atomic interaction)
  - Exchange correction (negligible for β⁺ but check)

### B5. Neutrino Spectrum Support

- [ ] Add convenience analyser methods for plotting neutrino spectra
  - Convert electron spectrum to neutrino spectrum
  - Plot both on same figure with proper scaling
- [ ] Adopt corrections that differ for neutrino spectrum
  - Implement neutrino radiative correction h(W̃, W0) per Sirlin 2011
  - Conversion formula: f_ν(W_ν) = f_e(W̃) · [1 + (α/2π)(h − g)]
  - Neutrino-specific corrections (no Fermi function, different radiative terms)

### B6. Neutrino Mass Effects

- [ ] Implement option to produce spectra with non-zero neutrino mass
  - Modify phase space: (W0 − W)² → (W0 − W)² · √(1 − m_ν²/(W0−W)²)
  - Add kink feature near endpoint for spectral distortion
  - Useful for neutrino mass search applications
- [ ] Research and document implications of significant neutrino mass
  - Which corrections need m_ν as input parameter?
  - Impact on radiative corrections (mass singularities)
  - Impact on endpoint fitting procedures

### B7. Automated Report Generation

- [ ] Implement extensive procedurally generated PDF report
  - Comprehensive details of the calculation
  - Component plots (phase space, Fermi, each correction factor)
  - Total spectrum with all corrections applied
  - Residual analysis, fit quality metrics
  - Parameter table with uncertainties
  - References to theoretical sources
- [ ] Extend report generation to fitters
  - Fit results with covariance matrix
  - Parameter correlations
  - Confidence intervals
  - Goodness-of-fit statistics

### B8. CLI & Output Improvements

- [x] Implement structured logging system
  - Log to stdout with configurable verbosity levels (-v=INFO, -vv=DEBUG)
  - Visualize workflow steps: data loading, parameter resolution, calculation progress
  - Log all used parameters: Z_parent, Z_daughter, A_number, endpoint, transition type, enabled corrections
  - Optional log file output with timestamped filename
  - Integrate Python `logging` module with custom formatter
- [x] Add metadata header to CSV output files
  - Date and time of calculation (ISO 8601 format)
  - Nuclide information: parent/daughter symbol, Z, A
  - Decay/transition info: endpoint energy, transition type, forbiddenness
  - Calculation parameters: energy step, enabled corrections, detector settings
  - Software version and git commit hash (7 chars)
  - Format: YAML-style comments at top of CSV file
- [ ] Optimize CLI argument design and add sanity checks
  - Auto-deduce `transition_type` from `decay_type` and nuclear data — remove as explicit CLI parameter
  - Auto-deduce `decay_type` from Z_parent vs Z_daughter difference — remove as explicit CLI parameter
  - Add sanity checks for custom input:
    - `|Z_parent - Z_daughter|` must equal 1 (beta decay)
    - For `beta_minus`: Z_daughter must be Z_parent + 1
    - For `beta_plus`/`ec`: Z_daughter must be Z_parent - 1
    - Warn if decay_type contradicts Z values (e.g., beta+ with Z_daughter > Z_parent)
    - Validate endpoint_MeV > 0 and endpoint_MeV > level_energy_keV/1000
    - Cross-check transition_type against ENSDF forbiddenness (if paceENSDF source)
- [x] Add `--dry-run` option to validate input and display resolved parameters without calculation
- [x] Add `--version` flag
- [x] Add `-q/--quiet` flag to suppress terminal output
- [x] Add `./output/` directory for debug verification artifacts (gitignored except `.gitkeep`)
- [x] Add debug verification step to development workflow: `bs_pnpi -vv` with parameter consistency checks across all components, documented in DEVELOPMENT.md §4.5
- [x] Add plot output with ID tracking (commit hash + UTC timestamp) and two display modes:
  - Normal mode (default): single spectrum plot with nuclear data header
  - Debug mode (-vv): 4-panel view with all correction factors and deviations

______________________________________________________________________

## Deferred (Not Planned)

### Additional Corrections

Not implementing at this time. Additional corrections (weak magnetism, recoil kinematics, etc.) require a complex implementation with project-wide consequences. This is a far-future plan, irrelevant to the current ⁹⁹Tc analysis goal where the shape factor C(W) absorbs all unmodelled energy-dependent effects.

### Nuclear Shape Factors C(Z,W)

Not implementing theoretical shape factor calculations. The entire point of this project is to **extract** C(W) from experimental data, not to compute it from nuclear theory (shell model, QRPA, etc.). Theoretical shape factor calculations require:

- Nuclear structure models (NUSHELLX, shell model codes)
- One-body transition densities (OBTDs)
- Coulomb displacement energy calculations
- Multiple effective interactions

These are handled by the nuclear physics community. Our role is to provide the analysis framework that extracts C(W) from data, then compare the extracted shape to theoretical predictions.

______________________________________________________________________

## Current Status

**Version:** 0.3.1 → **In progress: 0.4.0 — Package restructuring**

### Completed (v0.3.x)
- All theoretical corrections: phase space, Fermi function, finite size, screening, exchange, radiative (with delta_cut resummation)
- Detector response module with analytical models (Gaussian, Gaussian+tail, Tikhonov), convolution API
- χ² curve fitting framework (CurveFitter) with confidence intervals, profile likelihood, correlation analysis
- C(W) shape factor extraction pipeline (CWExtractor) with Kurie plot analysis, parametrized fitting, g_V/g_A extraction
- Multi-branch decay support (unified single/multi-branch)
- CLI interface (`bs_pnpi`) with paceENSDF integration (`--nuclide`), structured logging (-v/-vv/-q), --dry-run, --version, --log-file
- CSV metadata headers with branch info
- JSON input with full detector param support
- paceENSDF integration: `get_decay_info_from_paceENSDF()` returns DecayInfo with Z, A, endpoint, transition type, half-life, spin/parity, branches
- Comprehensive test suite (262 tests, including 6 paceENSDF integration tests)
- Notebook quality control with nbmake and auto-save plot hooks
- Plot output: unified multi-branch approach, externalized headers, compact legends, consistent titles

### In Progress (v0.4.0)
- **Track B: Package restructuring** — separating into `beta_spectrum` (theory), `exp_data` (experimental data), `fitter` (analysis framework)
  - Phase 1: Create `exp_data` module (raw data loading, calibration, simple fitters, corrections)
  - Phase 2: Create `fitter` module (theoretical model + convolution, multi-parameter fitting, C(W) extraction)
  - Phase 3: Clean up `beta_spectrum` (remove CWExtractor, CurveFitter, detector convolution)
  - Phase 4: Update package structure, imports, tests, documentation

### Next Immediate Step
- **v0.4.0 Phase 1:** Create `exp_data/` module with `raw_data.py`, `calibration.py`, `fitters.py`, `corrections.py`, `spectrum.py`
