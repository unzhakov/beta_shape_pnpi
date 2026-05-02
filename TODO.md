# Strategic Development Roadmap

## Scientific Goal

Extract the shape factor C(W) from the experimentally measured beta spectrum of ⁹⁹Tc, and parametrize C(W) via the vector coupling constant g_V from the theory of beta decay.

This requires repeating the analysis pipeline of [Paulsen et al., Phys. Rev. C **110**, 055503 (2024)](https://doi.org/10.1103/PhysRevC.110.055503), which used Metallic Magnetic Calorimeters (MMCs) to measure the ⁹⁹Tc beta spectrum with sub-keV threshold and ~100 eV resolution.

See `docs/refs/2024_Paulsen_T99_beta_spectrum/summary.md` for full analysis details.

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

- [ ] Implement analytical detector response function for 4π semiconductor detector
  - Gaussian core with low-energy tail (charge collection effects)
  - Parameters: energy resolution σ(E), tail fraction, tail shape
  - Energy-dependent resolution σ(E) = a + b·√E
- [ ] Implement convolution/deconvolution routines
  - Convolve theoretical spectrum with detector response for comparison to data
  - Deconvolve measured spectrum to recover true shape (for initial C(W) extraction)
  - Iterative unfolding (Richardson-Lucy or similar)
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

- [ ] Implement fitter routine to extract experimental C(W) from data
  - χ² minimization: theoretical spectrum (with all corrections) vs measured spectrum
  - Free parameters: C(W) shape parameters, endpoint energy, normalization, background
  - Covariance matrix and uncertainty propagation
- [ ] Parametrize C(W) in terms of g_V (and g_A)
  - Fit C(W) data to theoretical parametrization
  - Extract g_V^eff and g_A^eff from the fit
  - Compare with Paulsen et al. results: g_A^eff = 0.574(36), g_V^eff = 0.376(5)
- [ ] Systematic uncertainty analysis
  - Vary correction implementations within uncertainties
  - Test sensitivity to detector response model
  - Background model dependence

______________________________________________________________________

## Track B: Package Usability Improvements

### B1. Nuclear Data Integration

- [ ] Add usage of `paceENSDF` package to retrieve nuclear data for decays
  - Q-values, half-lives, branch intensities, decay schemes
  - Automated setup of W0, Z, A parameters from ENSDF
  - Support for multiple isotopes without manual input

### B2. Input Flexibility

- [ ] Generalize parameter inputs via custom input file (YAML/JSON)
  - Declare isotope, transition type, detector parameters
  - Toggle corrections on/off per-component
  - Reproducible calculation configurations
- [ ] Retrieve parameters directly from ENSDF database
  - Auto-populate W0, Z, A, transition type from isotope name
  - Reduce manual configuration errors

### B3. Multiple Branch Decays

- [ ] Add support for multiple branch decays
  - Calculate total spectral shape weighted by branch intensities
  - Handle decays with significant branching to excited states
  - Sum contributions from all branches with proper weighting

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

______________________________________________________________________

## Deferred (Not Planned)

### Nuclear Recoil Correction (`recoil.py`)

Not implementing at this time. The recoil correction (weak magnetism, recoil kinematics) is a small O(q/M_N) effect (~0.1%) that is subdominant to the shape factor C(W) we aim to extract. In the ⁹⁹Tc analysis, the shape factor absorbs all unmodelled energy-dependent effects — including recoil — so implementing recoil explicitly is unnecessary for the scientific goal.

### Nuclear Shape Factors C(Z,W)

Not implementing theoretical shape factor calculations. The entire point of this project is to **extract** C(W) from experimental data, not to compute it from nuclear theory (shell model, QRPA, etc.). Theoretical shape factor calculations require:

- Nuclear structure models (NUSHELLX, shell model codes)
- One-body transition densities (OBTDs)
- Coulomb displacement energy calculations
- Multiple effective interactions

These are handled by the nuclear physics community. Our role is to provide the analysis framework that extracts C(W) from data, then compare the extracted shape to theoretical predictions.

______________________________________________________________________

## Current Status

**Version:** 0.2.0\
**Implemented:** Phase space, Fermi function, finite size, screening, exchange, radiative corrections (with delta_cut resummation). Detector response module with analytical models (Gaussian, Gaussian+tail, Tikhonov), convolution API, declarative config integration. Comprehensive test suite (106 tests). Notebook quality control with nbmake and auto-save plot hooks.

**Completed:** A2 — detector response function and convolution routines.

**Next immediate step:** A4 — implement fitter routine to extract experimental C(W) from data.
