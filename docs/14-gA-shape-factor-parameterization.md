---
title: "14 – gA-Driven Shape Factor Parameterization"
date: 2026-05-13
tags:
  - beta-spectrum/shape-factor
  - beta-spectrum/gA
  - beta-spectrum/forbidden
status: active
aliases: [Spectrum Shape Method, SSM, gA quenching, shape factor parameterization]
related_components: [[10-nuclear-structure]], [[07-recoil-effects]]
---

# gA-Driven Shape Factor Parameterization — The Spectrum Shape Method (SSM)

## Overview

For **forbidden nonunique β decays**, the shape factor $C(Z,W)$ depends on the weak axial-vector coupling constant $g_A$ in a nontrivial, energy-dependent way. This creates a powerful observable: by fitting the theoretical spectrum shape to high-precision experimental data, one can **extract the effective value of $g_A$** (and $g_V$) for a given nuclear transition. This technique is called the **Spectrum Shape Method (SSM)**, introduced in Kostensalo et al. (2017).

> [!tip] Why this matters for 99Tc
> The $^{99}\text{Tc}$ β spectrum is a **second-forbidden nonunique** transition ($9/2_1^+ \to 5/2_1^+$). Its shape factor is **strongly sensitive** to $g_A$ — making it an excellent laboratory for determining the quenched effective coupling constants. Paulsen et al. (2024) measured the spectrum with MMCs down to $<1$ keV and extracted $g_V^{\text{eff}} = 0.376(5)$, $g_A^{\text{eff}} = 0.574(36)$.

## 1. Theoretical Foundation

### 1.1 Shape Factor Decomposition

The shape factor $C(w_e)$ decomposes into vector, axial-vector, and mixed parts:

$$C(w_e) = g_V^2\, C_V(w_e) + g_A^2\, C_A(w_e) + g_V g_A\, C_{VA}(w_e) \tag{1}$$

where $w_e = W$ is the total electron energy in units of $m_e c^2$. The three components $C_V$, $C_A$, $C_{VA}$ are **purely nuclear-structure dependent** — they are computed from nuclear matrix elements (NMEs) and are **independent of the coupling constants**.

After integration over the spectrum:

$$\tilde{C} = g_V^2\, \tilde{C}_V + g_A^2\, \tilde{C}_A + g_V g_A\, \tilde{C}_{VA} \tag{2}$$

The integrated components $\tilde{C}_i$ are energy-independent constants.

### 1.2 The β Spectrum Equation

The full differential spectrum (Paulsen et al., 2024, Eq. 5):

$$N(W)\,dW = \frac{G_F^2 \cos^2\theta_C}{2\pi^3}\, F(Z,W)\, p W (W_0 - W)^2\, X(W)\, C(W)\, \rho(Z,W)\, dW$$

where:
- $F(Z,W)$ — Fermi function (Coulomb correction)
- $X(W)$ — screening and exchange corrections
- $C(W)$ — **shape factor** (nuclear structure)
- $\rho(Z,W)$ — atomic overlap correction
- $G_F \cos\theta_C$ — effective Fermi constant

### 1.3 Traditional Parametrization (Ad-Hoc)

Historically, the shape factor for $^{99}\text{Tc}$ has been approximated as **first-forbidden unique** (Reich & Schüpferling, 1974; standard reference [15]):

$$C_{\text{adhoc}}(W) = q^2 + \lambda p^2 = (W_0 - W)^2 + \lambda (W^2 - 1) \tag{3}$$

with $\lambda$ as a free fitting parameter. This is **not derived from nuclear structure** — it is an empirical ansatz. Paulsen et al. (2024) extracted $\lambda \approx 0.65$ from fits, but this parametrization:
- Cannot distinguish $g_V$ from $g_A$ effects
- Provides no insight into nuclear structure
- Fails below $\sim$100 keV when compared to high-precision data

### 1.4 SSM Parametrization (Theory-Driven)

The SSM replaces Eq. (3) with the **microscopically computed** shape factor:

$$C_{\text{SSM}}(W; g_A) = g_V^2\, C_V(W) + g_A^2\, C_A(W) + g_V g_A\, C_{VA}(W) \tag{4}$$

where $C_V(W)$, $C_A(W)$, $C_{VA}(W)$ are computed from shell-model (or MQPM) one-body transition densities (OBTDs).

## 2. Algorithm for Parameterization

### 2.1 Input Requirements

| Quantity | Source | Notes |
|---|---|---|
| $C_V(W)$, $C_A(W)$, $C_{VA}(W)$ | Nuclear shell model (NuShellX@MSU) or MQPM | Computed for the specific transition |
| $W_0$ (endpoint) | Independent measurement or fit | Paulsen et al.: $W_0 = 295.82(16)$ keV |
| Experimental spectrum $N_{\text{exp}}(W)$ | Calorimetry / semiconductor / MMC | Corrected for detector response, background, energy losses |
| $F(Z,W)$, $X(W)$, $\rho(Z,W)$ | Calculator modules | Already implemented in project |
| $g_V$ | Fixed (CVC: $g_V = 1.0$, or treated as free) | Paulsen et al. used free $g_V$ fit |

### 2.2 Computation of Shape Factor Components

The shape factor in the **Behrens-Bühring** formalism (truncated lepton current approximation, Paulsen et al. Eq. 10):

$$C(W) = \sum_{K,k_e,k_\nu} \left[ e_k^2 W k_e\, m_K(k_e,k_\nu) - 2\mu e_k W k_e\, \gamma\, m_K(k_e,k_\nu) \right]$$

where the sum runs over lepton quantum numbers $k_e$, $k_\nu$ and nuclear multipole order $K$. The key point: **each component** $C_i(W)$ ($i = V, A, VA$) is computed separately from the NMEs.

**Practical computation via shell model:**

1. **OBTD calculation** (NuShellX):
   - Compute one-body transition densities for the transition
   - For $^{99}\text{Tc}$: dominant transition $2d_{5/2} \to 1g_{9/2}$ with OBTD $\approx 0.478$
   - Sub-dominant: $1g_{7/2} \to 1g_{9/2}$ with OBTD $\approx 0.010$

2. **Coulomb displacement energy** $E_C$ (for relativistic vector current):
   - Leading order: $E_C^{(3)} = \frac{\int_0^\infty g_f V(r) g_i (r/R)^K r^2 dr}{\int_0^\infty g_f g_i (r/R)^K r^2 dr}$
   - Energy-dependent: $E_C^{(4)}$ includes full lepton current dependence
   - For $^{99}\text{Tc}$: $E_C^{(4)} \approx 10.5$ MeV (GLEKPN valence space)

3. **Vector current from CVC**:
   $$V_F^{221} = -\sqrt{\frac{R}{10}}\left[W_0 - (m_n - m_p) + E_C\right] V_F^{220}$$
   For $^{99}\text{Tc}$, $[W_0 - (m_n - m_p)] = -0.487$ MeV — small Q-value makes $E_C$ critical.

### 2.3 Fitting Procedure

```
Step 1: Precompute shape factor components C_i(W) from nuclear model
        → Store as lookup tables: C_V(W), C_A(W), C_VA(W) for W ∈ [1, W0]

Step 2: For trial (gA, gV) pair:
        a. Compute C(W; gA, gV) = gV²·C_V(W) + gA²·C_A(W) + gV·gA·C_VA(W)
        b. Compute full spectrum:
           N_theo(W) ∝ p·W·(W0−W)²·F(Z,W)·X(W)·C(W; gA, gV)·ρ(Z,W)
        c. Normalize: ∫ N_theo(W) dW = 1

Step 3: Compare with experimental spectrum:
        χ²(gA, gV) = Σ [N_exp(W_i) − N_theo(W_i; gA, gV)]² / σ_i²

Step 4: Minimize χ² over (gA, gV) grid
        → Extract best-fit values and confidence intervals
```

### 2.4 Energy Range Selection

Paulsen et al. (2024) recommend fitting in the range **20–275 keV** to avoid:
- **Below ~20 keV**: detector threshold effects, atomic exchange modeling uncertainties
- **Above ~275 keV**: endpoint region with poor statistics, resolution distortion

The χ² landscape is approximately quadratic, enabling simple uncertainty estimation at χ² = χ²_min + 1.

## 3. Results from Literature

### 3.1 $^{99}\text{Tc}$ — Paulsen et al. (2024)

From the MMC measurement (PTB/LNHB), the extracted values are:

| Quantity | Value | Method |
|---|---|---|
| $Q_\beta$ | 295.82(16) keV | Fit to spectrum shape |
| $g_V^{\text{eff}}$ | 0.376(5) | From half-life normalization |
| $g_A^{\text{eff}}$ | 0.574(36) | From spectrum shape fit |
| $\bar{E}_\beta$ | 98.51(23) keV | From theoretical spectrum |
| $\log f$ | −0.476 60(22) | From phase space integration |
| $\log ft$ | 12.3478(23) | Combined |

**Key finding**: The shape-fit alone gives $g_A^{\text{eff}} = 1.526(92)$ (enhanced, not quenched), but when the half-life is used to normalize the overall rate, the **quenched** values $g_V^{\text{eff}} = 0.376$, $g_A^{\text{eff}} = 0.574$ emerge — consistent with first-forbidden nonunique transitions (Suhonen review: $g_A^{\text{eff}} \approx 0.9$, $g_V^{\text{eff}} \approx 0.3$–0.7).

### 3.2 $^{87}\text{Rb}$ — Kostensalo et al. (2017)

Third-forbidden nonunique ($3/2^- \to 9/2^+$):
- Strong $g_A$ dependence in spectrum shape
- MQPM and NSM agree on shape evolution
- $g_A \approx 0.9$ gives best match to available data
- Best candidate for $g_A$ determination in third-forbidden decays

### 3.3 $^{113}\text{Cd}$ and $^{115}\text{In}$ — Kostensalo et al. (2017, 2017a)

Fourth-forbidden nonunique:
- Most dramatic $g_A$ dependence: bell-shaped spectrum at $g_A \approx g_V$, monotonic otherwise
- Consistent $g_A \approx 0.9$ extracted across three nuclear models (MQPM, NSM, IBFM-2)
- Demonstrates **robustness of SSM** against nuclear model details

### 3.4 Summary of $g_A$ Sensitivity by Forbiddenness

| Transition | Forbiddenness | $g_A$ Sensitivity | SSM Viability |
|---|---|---|---|
| $^{99}\text{Tc}$ | 2nd nonunique | **Strong** | ✅ Excellent |
| $^{94}\text{Nb}$, $^{98}\text{Tc}$ | 2nd nonunique | **Strong** | ✅ Good |
| $^{87}\text{Rb}$ | 3rd nonunique | **Strong** | ✅ Good |
| $^{113}\text{Cd}$, $^{115}\text{In}$ | 4th nonunique | **Very strong** | ✅ Excellent |
| $^{137}\text{Cs}$ (to $11/2^-$) | 1st unique | Weak (NLO only) | ⚠️ Limited |
| $^{137}\text{Cs}$ (to $3/2^+$) | 2nd nonunique | **None** (cancellation) | ❌ Not applicable |
| First-forbidden nonunique | 1st nonunique | **Weak** | ❌ Not useful |

**Systematic pattern**: Even-forbidden decays show strong $g_A$ dependence; odd-forbidden decays are mostly insensitive.

## 4. Implementation Plan

### 4.1 New Module: `components/shape_factor_gA.py`

```
ShapeFactorGA
├── __init__: load C_V(W), C_A(W), C_VA(W) from precomputed tables or nuclear model
├── evaluate(W, gA, gV=1.0): compute C(W; gA, gV)
├── evaluate_components(W): return (C_V, C_A, C_VA) separately
├── chi2_spectrum(gA, gV, W_grid, N_exp, sigma): χ² for spectrum fit
├── fit_gA(W_grid, N_exp, sigma, W_min, W_max): minimize χ², return (gA_best, gA_err)
├── fit_gA_gV(W_grid, N_exp, sigma, W_min, W_max): fit both, return (gA, gV, cov)
└── mean_energy(gA, gV, W0): compute ⟨E⟩ from theoretical spectrum
```

### 4.2 Integration with BetaSpectrum

The shape factor becomes a new `SpectrumComponent`:

```python
class ShapeFactorGA(SpectrumComponent):
    """gA-dependent shape factor for forbidden nonunique decays."""
    
    def __init__(self, C_V_table, C_A_table, C_VA_table, gA=1.27, gV=1.0):
        # Load precomputed shape factor components
        # C_i tables: W array → C_i(W) values
    
    def evaluate(self, W):
        return gV**2 * self.C_V(W) + gA**2 * self.C_A(W) + gV*gA * self.C_VA(W)
    
    def integrate(self, W0):
        """Compute integrated shape factor for half-life prediction."""
        return integrate(C(W) * p * W * (W0-W)**2 * F(Z,W) dW, 1, W0)
```

### 4.3 Precomputed Tables from Nuclear Model

For $^{99}\text{Tc}$, the shape factor components should be precomputed from NuShellX@MSU (GLEKPN valence space, Mach interaction):

| W (keV) | $C_V(W)$ | $C_A(W)$ | $C_{VA}(W)$ |
|---|---|---|---|
| 1.001 | ... | ... | ... |
| ... | ... | ... | ... |
| 1.581 | ... | ... | ... |

(Endpoint: $W_0 = 1 + 295.82/510.998950 \approx 1.579$)

### 4.4 Fitting Pipeline

```python
def fit_gA_for_Tc95(data_file, W_min=20*keV, W_max=275*keV):
    """Fit gA to 99Tc beta spectrum data."""
    W, N_exp, sigma = load_spectrum(data_file)
    
    # Load precomputed shape factor components
    sf = ShapeFactorGA("Tc95_CV_CA_CVA_tables.npy")
    
    # Grid search + refinement
    best_chi2 = inf
    for gA in np.linspace(0.4, 1.8, 141):
        chi2 = sf.chi2_spectrum(gA=gA, gV=1.0, W=W, N_exp=N_exp, sigma=sigma,
                                W_min=W_min, W_max=W_max)
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_gA = gA
    
    # Refine around minimum
    gA_opt = optimize.minimize_scalar(...)
    
    # Uncertainty from χ² = χ²_min + 1
    gA_low, gA_high = find_chi2_one_boundary(gA_opt, sf, ...)
    
    return gA_opt, gA_low, gA_high
```

## 5. Uncertainty Budget (from Paulsen et al., 2024)

| Component | Uncertainty on $g_A^{\text{eff}}$ | Relative |
|---|---|---|
| Fit method (statistics) | 0.0745 | 4.9% |
| Coulomb displacement energy | 0.0425 | 2.8% |
| Nuclear model | 0.0294 | 1.9% |
| Lepton current treatment | 0.0127 | 0.8% |
| Fit range | 0.0061 | 0.4% |
| Maximum energy | 0.0038 | 0.2% |
| Radiative corrections | 0.0033 | 0.2% |
| Atomic exchange | 0.0029 | 0.2% |
| **Combined** | **0.0919** | **6.0%** |

Dominant uncertainty: **fit method / statistics**. The nuclear model uncertainty is surprisingly small (1.9%), confirming SSM robustness.

## 6. Key References

| Ref | Citation | Relevance |
|---|---|---|
| [SSM-1] | Kostensalo, Haaranen, Suhonen, *Phys. Rev. C* **95**, 044313 (2017) | Original SSM paper; 26 decays studied |
| [SSM-2] | Kostensalo, Suhonen, *Phys. Rev. C* **96**, 024317 (2017) | gA-driven evolution; 16 decays; Table I with integrated shape factors |
| [SSM-3] | Paulsen et al., *Phys. Rev. C* **110**, 055503 (2024) | $^{99}\text{Tc}$ high-precision MMC measurement; $g_A$ extraction |
| [SSM-4] | Kostensalo et al., *Phys. Rev. C* **95**, 024327 (2017) | $^{113}\text{Cd}$ analysis with three nuclear models; SSM robustness |
| [Hayen] | Hayen et al., *Rev. Mod. Phys.* **90**, 015008 (2017) | Comprehensive review; BB formalism; Eq. 79–81 |
| [Suhonen] | Suhonen, *From Nucleons to Nucleus* (Springer, 2007) | Nuclear structure background; effective coupling constants |
| [Behrens] | Behrens & Bühring, *Electron Radial Wave Functions and Nuclear Beta Decay* (1982) | BB formalism reference; tabulated $I(k_e,m,n,\rho)$ functions |

## 7. Connection to Existing Modules

This module integrates with:
- **`components/fermi.py`** — Fermi function $F(Z,W)$
- **`components/exchange.py`** — Exchange correction $X(W)$ (extended to forbidden transitions)
- **`components/finite_size.py`** — Finite-size correction (daughter nucleus radius)
- **`components/screening.py`** — Screening correction
- **`spectrum.py`** — `BetaSpectrum` multiplicative composition
- **`docs/10-nuclear-structure.md`** — Theoretical background on shape factors

## 8. Practical Considerations

### 8.1 Which transitions to prioritize?

For $g_A$ extraction via SSM, the selection criteria are:
1. **Strong $g_A$ dependence** in spectrum shape (even-forbidden, nonunique)
2. **High branching ratio** (preferably > 90%)
3. **Measurable Q-value** with good energy resolution
4. **Low background** environment

For the $^{99}\text{Tc}$ project, this is the primary candidate — it satisfies all criteria.

### 8.2 Quenching vs. Enhancement

The extracted $g_A^{\text{eff}}$ from shape fits alone can be **enhanced** ($>1.27$) because the shape factor normalization absorbs the overall rate. The **quenched** values emerge only when the half-life is used for absolute normalization:

$$g_i^{\text{eff, quenched}} = g_i^{\text{eff, shape}} \times \sqrt{\frac{t_{1/2}^{\text{theo}}}{t_{1/2}^{\text{exp}}}}$$

For $^{99}\text{Tc}$: $t_{1/2}^{\text{theo}} \approx 30 \times 10^3$ a vs. $t_{1/2}^{\text{exp}} \approx 212 \times 10^3$ a → quenching factor $\approx 0.38$.

### 8.3 Nuclear Model Dependencies

While SSM is robust against nuclear model details (MQPM vs. NSM vs. IBFM-2 agree to $\sim$10%), the **absolute normalization** of shape factor components varies:
- For $^{99}\text{Tc}$: MQPM gives $\sim$2× larger $\tilde{C}$ than NSM
- This affects half-life predictions but **not** the shape evolution
- The $g_A$ extraction from shape is thus model-independent to $\sim$2% (see uncertainty budget)

## 9. Future Extensions

1. **First-forbidden unique decays**: Include NLO terms for $g_A$ sensitivity
2. **Highly forbidden decays** (4th–6th order): $^{96}\text{Zr}$, $^{48}\text{Ca}$ — strong $g_A$ dependence in even-forbidden channels
3. **Two-body currents**: Recent ab initio calculations suggest quenching may be unnecessary when 2BC are included (Carbone et al.)
4. **Multi-parameter fits**: Simultaneous fit of $g_A$, $g_V$, $Q_\beta$, and nuclear structure parameters
5. **Machine learning surrogate**: Train emulator on shape factor components for rapid χ² evaluation
