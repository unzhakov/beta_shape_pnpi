---
title: "SSM — Spectrum Shape Method for gA Extraction"
date: 2026-05-14
tags:
  - gA/extraction/algorithm
  - gA/ssm
  - beta-spectrum/forbidden
status: active
aliases: [Spectrum Shape Method, gA fitting, shape factor fit]
related_notes:
  - gA-extraction/02-ssnme
  - gA-extraction/03-tc99-case-study
  - gA-extraction/04-implementation-todo
  - nuclear-structure/01-shape-factors
  - corrections/05-radiative
---

# Spectrum Shape Method (SSM) — Step-by-Step Algorithm for gA Extraction

## 1. Physical Basis

For **forbidden nonunique β decays**, the shape factor $C(w_e)$ decomposes into three energy-dependent components with distinct coupling-constant dependence:

$$C(w_e) = g_V^2\, C_V(w_e) + g_A^2\, C_A(w_e) + g_V g_A\, C_{VA}(w_e) \tag{1}$$

where $w_e$ is the total electron energy in units of $m_e c^2$. The components $C_V$, $C_A$, $C_{VA}$ are computed from nuclear many-body theory and are **independent of the coupling constants**. The spectrum shape is sensitive to $g_A$ through the interference between vector, axial-vector, and mixed terms — a purely nuclear-structure effect unique to nonunique forbidden decays.

> [!note] Why nonunique decays?
> Unique forbidden decays have a single dominant NME → universal spectrum shape independent of nuclear structure. Nonunique decays involve multiple NMEs → the relative weight of vector vs. axial contributions depends on $g_A$, creating observable spectral shape differences.

## 2. The Full Spectrum Model

The differential spectrum (Paulsen et al., 2024, Eq. 5):

$$N(W) = \mathcal{N} \cdot F(Z,W) \cdot p W (W_0 - W)^2 \cdot X(W) \cdot S(Z,W) \cdot R(W) \cdot C(W; g_A, g_V) \cdot \rho(Z,W)$$

where:
- $\mathcal{N}$ — overall normalization (absorbs $G_F^2 \cos^2\theta_C / 2\pi^3$)
- $F(Z,W)$ — Fermi function (Coulomb correction) — [[corrections/01-fermi-function]]
- $p W (W_0 - W)^2$ — phase space — [[02-phase-space]]
- $X(W)$ — atomic exchange correction — [[corrections/04-exchange]]
- $S(Z,W)$ — atomic screening — [[corrections/03-screening]]
- $R(W)$ — outer radiative correction — [[corrections/05-radiative]]
- $C(W; g_A, g_V)$ — **shape factor** — [[nuclear-structure/01-shape-factors]]
- $\rho(Z,W)$ — atomic overlap correction — [[corrections/08-atomic-overlap]]

## 3. The Algorithm

### Step 1: Compute shape factor components from nuclear model

**Input:** Nuclear many-body model (NSM, MQPM, etc.) with appropriate valence space and Hamiltonian for the transition of interest.

**Output:** Three energy-dependent lookup tables: $C_V(W)$, $C_A(W)$, $C_{VA}(W)$ for $W \in [1, W_0]$.

**Method:**
1. Run shell model (or MQPM) to obtain one-body transition densities (OBTDs)
2. For each OBTD, compute the contribution to each shape factor component using the Behrens-Bühring formalism (truncated lepton current approximation)
3. Sum over all contributing multipole orders $K$ and lepton quantum numbers $k_e, k_\nu$

**CVC constraint for vector current:**
The relativistic vector current includes the small relativistic NME (sNME), which is related to the large vector NME (l-NME) via the CVC hypothesis:

$$\text{sNME} = \frac{1}{2} \left[ W_0 - (m_n - m_p) + E_C \right] \cdot \text{l-NME}$$

where $E_C$ is the Coulomb displacement energy. The sNME is critical for small-$Q$ decays like $^{99}\text{Tc}$ ($[W_0 - (m_n - m_p)] = -0.487$ MeV).

**Practical note:** For $^{99}\text{Tc}$, the dominant transition is $2d_{5/2} \to 1g_{9/2}$ with OBTD $\approx 0.478$. The sub-dominant $1g_{7/2} \to 1g_{9/2}$ transition contributes OBTD $\approx 0.010$. See Kostensalo & Suhonen (2017), Table VII.

### Step 2: Build theoretical spectrum for trial $g_A$

For each trial value of $g_A$ (with $g_V = 1.0$ fixed by CVC):

1. **Compute shape factor:**
   $$C(W; g_A) = C_V(W) + g_A^2\, C_A(W) + g_A\, C_{VA}(W)$$
   (using $g_V = 1.0$)

2. **Compute full spectrum:**
   $$N_{\text{theo}}(W; g_A) = \mathcal{N} \cdot F(Z,W) \cdot p W (W_0 - W)^2 \cdot X(W) \cdot S(Z,W) \cdot R(W) \cdot C(W; g_A) \cdot \rho(Z,W)$$

3. **Normalize:** Choose $\mathcal{N}$ so that $\int_1^{W_0} N_{\text{theo}}(W; g_A)\, dW = 1$ (or match the total counts in the experimental spectrum)

### Step 3: Compare with experimental spectrum

**Input:** Experimental spectrum data $(W_i, N_i, \sigma_i)$ for $i = 1, \dots, N_{\text{bins}}$.

**χ² function:**
$$\chi^2(g_A) = \sum_{i \in \text{fit range}} \frac{\left[N_i - N_{\text{theo}}(W_i; g_A)\right]^2}{\sigma_i^2}$$

**Fit range selection (Paulsen et al., 2024):**
- **Lower bound:** ~20 keV — below this, detector threshold effects and atomic exchange modeling dominate
- **Upper bound:** ~275 keV for $^{99}\text{Tc}$ — above this, endpoint region with poor statistics and resolution distortion
- Rationale: the $g_A$-driven shape differences are most pronounced in the mid-energy range

### Step 4: Minimize χ² and extract $g_A$

**Grid search:** Evaluate $\chi^2(g_A)$ on a grid (e.g., $g_A \in [0.4, 1.8]$ with step 0.01).

**Refinement:** Use `scipy.optimize.minimize_scalar` (Brent's method) around the grid minimum.

**Uncertainty:** Find $g_A$ values where $\chi^2 = \chi^2_{\text{min}} + 1$. These define the 68% confidence interval.

**Output:** $(g_A^{\text{best}}, g_A^{\text{low}}, g_A^{\text{high}})$

### Step 5 (optional): Simultaneous $g_A$ and $g_V$ fit

If $g_V$ is not fixed to 1.0 (e.g., testing CVC violation):

$$\chi^2(g_A, g_V) = \sum_i \frac{\left[N_i - N_{\text{theo}}(W_i; g_A, g_V)\right]^2}{\sigma_i^2}$$

Fit both parameters simultaneously. The correlation between $g_A$ and $g_V$ is typically strong.

### Step 6 (optional): Half-life normalization — quenched $g_A$

The $g_A$ extracted from **shape alone** can be enhanced ($>1.27$) because the normalization constant $\mathcal{N}$ absorbs the overall rate. To obtain the **quenched** values consistent with half-life:

1. Compute the theoretical partial half-life for the best-fit $(g_A, g_V)$:
   $$t_{1/2}^{\text{theo}} = \frac{\kappa}{\tilde{C}}, \quad \tilde{C} = \int_1^{W_0} C(W; g_A, g_V) \cdot p W (W_0 - W)^2 \cdot F(Z,W) \cdot X(W) \cdot S(Z,W) \cdot R(W) \, dW$$

2. Compare with experimental partial half-life:
   $$t_{1/2}^{\text{exp}} = \frac{t_{1/2}^{\text{total}}}{\text{branching ratio}}$$

3. Extract quenched values:
   $$g_i^{\text{quenched}} = g_i^{\text{shape}} \times \sqrt{\frac{t_{1/2}^{\text{theo}}}{t_{1/2}^{\text{exp}}}}$$

For $^{99}\text{Tc}$: $t_{1/2}^{\text{theo}} \approx 30 \times 10^3$ a vs. $t_{1/2}^{\text{exp}} \approx 212 \times 10^3$ a → quenching factor $\approx 0.38$.

## 4. The sNME Problem and ESSM

### The Problem

The sNME (small relativistic vector NME) is extremely difficult to compute in nuclear shell model because it involves contributions outside the valence shell. The CVC prediction gives a "canonical" value, but this assumes a perfect many-body theory — which finite valence spaces are not.

### The Solution: Enhanced SSM (ESSM)

Kumar et al. (2021) proposed treating the sNME as a **fitting parameter** alongside $g_A$. This is the **Enhanced Spectrum Shape Method (ESSM)**:

1. Fix $g_A$ at a trial value
2. Vary the sNME to simultaneously fit the **experimental half-life** (not just the shape)
3. Because the half-life depends quadratically on sNME, there are typically **two solutions** for each $g_A$ value:
   - **sNME(c):** closer to the CVC value
   - **sNME(f):** further from the CVC value
4. Compare both solutions' spectral shapes with the data

The "optimal" choice (closest to CVC) is theoretically motivated, but experimental validation is needed. See [[gA-extraction/02-ssnme]] for details.

## 5. Uncertainty Budget

| Component | Uncertainty on $g_A$ | Relative | Dominant for |
|---|---|---|---|
| Fit method (statistics) | 0.0745 | 4.9% | All SSM analyses |
| Coulomb displacement energy | 0.0425 | 2.8% | Small-$Q$ decays |
| Nuclear model | 0.0294 | 1.9% | Model comparison |
| Lepton current treatment | 0.0127 | 0.8% | High-$Z$ nuclei |
| Fit range | 0.0061 | 0.4% | All |
| Endpoint energy | 0.0038 | 0.2% | Precision $Q$-value needed |
| Radiative corrections | 0.0033 | 0.2% | All |
| Atomic exchange | 0.0029 | 0.2% | Low-energy region |
| **Combined** | **0.0919** | **6.0%** | — |

Source: Paulsen et al. (2024), Table VIII.

## 6. Key References

| Ref | Citation | Use |
|---|---|---|
| [SSM-orig] | Kostensalo, Haaranen, Suhonen, *Phys. Rev. C* **95**, 044313 (2017) | SSM introduction; 26 decays |
| [SSM-Tc94-98] | Kostensalo, Suhonen, *Phys. Rev. C* **96**, 024317 (2017) | gA-driven evolution; $^{94,98}$Tc, $^{99}$Tc |
| [Paulsen2024] | Paulsen et al., *Phys. Rev. C* **110**, 055503 (2024) | $^{99}$Tc MMC measurement; $g_A$ extraction (enhanced value) |
| [Song2025] | Song et al., arXiv:2506.17544 (2025) | $^{99}$Tc MMC measurement; quenched $g_A$ via ESSM |
| [ESSM] | Kumar, Srivastava, Suhonen, *Eur. Phys. J. A* **57**, 225 (2021) | Enhanced SSM formalism |
| [Ramalho2024] | Ramalho, Suhonen, *Front. Phys.* **12**, 1455778 (2024) | sNME analysis; five second-forbidden decays |
| [Hayen2017] | Hayen et al., *Rev. Mod. Phys.* **90**, 015008 (2017) | BB formalism; shape factor theory |
