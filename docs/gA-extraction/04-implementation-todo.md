---
title: "Implementation TODO — gA Extraction Pipeline"
date: 2026-05-14
tags:
  - gA/extraction/todo
  - implementation-plan
status: active
aliases: [gA pipeline, SSM implementation, shape factor module]
related_notes:
  - gA-extraction/01-ssm-algorithm
  - gA-extraction/02-ssnme
  - gA-extraction/03-tc99-case-study
---

# Implementation TODO — gA Extraction Pipeline

> [!info] Status
> **Research phase.** Experimental data to be added by user. Implementation to follow after documentation is complete.

## Phase 1: Nuclear Model Input (P0 — prerequisite)

- [ ] **1.1** Obtain precomputed shape factor components $C_V(W)$, $C_A(W)$, $C_{VA}(W)$ for $^{99}\text{Tc}$
  - [ ] Option A: Use Kostensalo & Suhonen (2017) MQPM results (published in PRC 96, 024317)
  - [ ] Option B: Use Ramalho & Suhonen (2024) NSM results (Front. Phys. 12, 1455778)
  - [ ] Option C: Run shell model with KShell or another available code
  - [ ] **Decision:** Which source to use? (See [[gA-extraction/03-tc99-case-study]] for comparison)
- [ ] **1.2** Store components as lookup tables (CSV or NumPy):
  - `data/Tc95_CV.csv` — $C_V(W)$ vs. $W$
  - `data/Tc95_CA.csv` — $C_A(W)$ vs. $W$
  - `data/Tc95_CVA.csv` — $C_{VA}(W)$ vs. $W$
  - Columns: `W (unitless)`, `C_V`, `C_A`, `C_VA`
  - $W$ range: $[1.001, 1.580]$ (corresponds to $T \in [0.5, 295]$ keV)
- [ ] **1.3** Decide on sNME handling:
  - [ ] Fix to CVC value (Paulsen et al. approach)
  - [ ] Fit as free parameter (ESSM — Song et al. approach) **RECOMMENDED**
  - If ESSM: also need l-NME from nuclear model

## Phase 2: Shape Factor Module (P1 — depends on 1.1)

Create `beta_spectrum/components/shape_factor_gA.py`:

- [ ] **2.1** `ShapeFactorGA` class:
  ```python
  class ShapeFactorGA:
      def __init__(self, C_V_table, C_A_table, C_VA_table, sNME=None):
          # Load tables, interpolate to arbitrary W
          # sNME: if None, use CVC value; if float, use fitted value
      
      def evaluate(self, W, gA, gV=1.0):
          """Compute C(W; gA, gV) = gV²·CV(W) + gA²·CA(W) + gV·gA·CVA(W)"""
      
      def evaluate_components(self, W):
          """Return (C_V(W), C_A(W), C_VA(W))"""
  ```
- [ ] **2.2** Interpolation: Use `scipy.interpolate.interp1d` with `kind='cubic'` for smooth evaluation
- [ ] **2.3** Validation: Check that $C(W)$ is monotonically reasonable for $g_A \in [0.4, 1.8]$
- [ ] **2.4** Unit test: `test_shape_factor_gA_consistency` — verify decomposition matches Eq. (1)

## Phase 3: Fitting Pipeline (P2 — depends on 1.2 + 2.1)

Create `beta_spectrum/gA_fit.py`:

- [ ] **3.1** `fit_gA_spectrum(W_exp, N_exp, sigma_exp, shape_factor, W_min, W_max)`:
  ```python
  def fit_gA_spectrum(W_exp, N_exp, sigma_exp, sf, W_min=0.020, W_max=0.275):
      """
      Fit gA to experimental spectrum using χ² minimization.
      
      Parameters
      ----------
      W_exp : array_like
          Electron total energy (unitless, W = 1 + T/m_ec²)
      N_exp : array_like
          Measured counts (or normalized spectrum)
      sigma_exp : array_like
          Uncertainties on N_exp
      sf : ShapeFactorGA
          Precomputed shape factor
      W_min, W_max : float
          Fit range in unitless W
      
      Returns
      -------
      gA_best, gA_low, gA_high, chi2_min, chi2_curve
      """
  ```
- [ ] **3.2** Grid search: $g_A \in [0.4, 1.8]$ with 0.01 step
- [ ] **3.3** Brent's method refinement around grid minimum
- [ ] **3.4** Confidence interval: find $g_A$ where $\chi^2 = \chi^2_{\text{min}} + 1$
- [ ] **3.5** χ² curve output for visualization
- [ ] **3.6** Unit test: `test_fit_gA_reproduces_input` — inject known $g_A$, verify recovery

## Phase 4: ESSM Extension (P2 — optional but recommended)

- [ ] **4.1** `fit_gA_sNME(W_exp, N_exp, sigma_exp, sf, l_NME, t_half_exp)`:
  - Simultaneous fit of $g_A$ and sNME
  - Match both spectral shape AND experimental half-life
  - Return both (sNME_c, sNME_f) solutions
- [ ] **4.2** Half-life computation:
  ```python
  def compute_half_life(sf, gA, gV, sNME, W0):
      """Compute t_1/2 from integrated shape factor"""
      C_integral = integrate(C(W; gA, gV, sNME) * phase_space * Fermi * ... dW)
      return kappa / C_integral
  ```
- [ ] **4.3** Validation against Song et al. (2025) results for $^{99}\text{Tc}$

## Phase 5: Validation (P3 — depends on 3.1 + 4.1)

- [ ] **5.1** Reproduce Paulsen et al. (2024) result: $g_A^{\text{eff}} = 1.526(92)$ with CVC-fixed sNME
- [ ] **5.2** Reproduce Song et al. (2025) result: $g_A^{\text{eff}} \approx 1.0\text{--}1.2$ with ESSM
- [ ] **5.3** Cross-check: polynomial fit $A \cdot W + B \cdot W^{-1} + C \cdot W^2$ vs. SSM for same $g_A$
- [ ] **5.4** Sensitivity analysis: vary $W_{\text{min}}$, $W_{\text{max}}$, nuclear model, sNME

## Phase 6: Experimental Data Analysis (P4 — depends on user providing data)

- [ ] **6.1** Import raw spectrum data (format TBD)
- [ ] **6.2** Apply detector response correction (if needed) — [[corrections/07-detector-response]]
- [ ] **6.3** Run `fit_gA_spectrum` with appropriate fit range
- [ ] **6.4** Run `fit_gA_sNME` (ESSM) if nuclear model input available
- [ ] **6.5** Generate publication-quality plots:
  - Measured spectrum + best-fit theoretical curve
  - Residuals plot
  - χ² vs. $g_A$ curve
  - Comparison with polynomial fit
- [ ] **6.6** Uncertainty budget: propagate detector, model, and systematic uncertainties

## Phase 7: Paper Draft Support

- [ ] **7.1** Generate tables: fitted $g_A$, sNME, $\chi^2$, confidence intervals
- [ ] **7.2** Generate figures for paper
- [ ] **7.3** Cross-reference with literature values
- [ ] **7.4** Draft methodology section using [[gA-extraction/01-ssm-algorithm]] as reference

## File Structure (target)

```
beta_shape_pnpi/
├── beta_spectrum/
│   ├── components/
│   │   └── shape_factor_gA.py        ← Phase 2
│   └── gA_fit.py                      ← Phase 3 + 4
├── data/
│   ├── Tc95_CV.csv                    ← Phase 1
│   ├── Tc95_CA.csv
│   ├── Tc95_CVA.csv
│   └── experimental/                  ← Phase 6 (user-provided)
├── tests/
│   └── test_shape_factor_gA.py        ← Phase 2 + 3
│   └── test_gA_fit.py
└── docs/
    └── gA-extraction/
        ├── 01-ssm-algorithm.md
        ├── 02-ssnme.md
        ├── 03-tc99-case-study.md
        └── 04-implementation-todo.md
```

## Decision Points

1. **Nuclear model source:** Which precomputed $C_i(W)$ to use? (MQPM vs. NSM)
2. **sNME strategy:** CVC-fixed or ESSM? (ESSM recommended)
3. **Fit range:** What $W_{\text{min}}$, $W_{\text{max}}$ for your detector?
4. **Detector corrections:** Do you need response deconvolution before fitting?

## References for Implementation

| Item | Source |
|---|---|
| Shape factor formalism | Kostensalo & Suhonen, PRC 96, 024317 (2017) |
| ESSM formalism | Kumar et al., EPJA 57, 225 (2021) |
| sNME analysis | Ramalho & Suhonen, Front. Phys. 12, 1455778 (2024) |
| $^{99}\text{Tc}$ MMC data | Paulsen et al., PRC 110, 055503 (2024) |
| $^{99}\text{Tc}$ ESSM analysis | Song et al., arXiv:2506.17544 (2025) |
| BB formalism | Behrens & Bühring (1982) |
