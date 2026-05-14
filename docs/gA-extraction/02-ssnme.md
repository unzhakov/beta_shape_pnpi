---
title: "sNME — Small Relativistic Vector NME"
date: 2026-05-14
tags:
  - gA/extraction/ssnme
  - beta-spectrum/forbidden
  - nuclear-structure/NME
status: active
aliases: [small relativistic vector nuclear matrix element, sNME, CVC vector current]
related_notes:
  - gA-extraction/01-ssm-algorithm
  - gA-extraction/03-tc99-case-study
  - nuclear-structure/01-shape-factors
---

# The Small Relativistic Vector NME (sNME) — Why It Matters

## What is the sNME?

In forbidden β decay, the relativistic vector current has two components:
- **Large NME (l-NME):** Computed from the large components of nucleon wave functions. Reliably calculable in NSM.
- **Small NME (sNME):** Couples the large and small components. **Extremely difficult** to compute because it involves contributions outside the valence shell.

In a perfect many-body theory with infinite valence space, the sNME is related to the l-NME by the **CVC hypothesis**:

$$\text{sNME} = \frac{1}{2} \left[ W_0 - (m_n - m_p) + E_C \right] \cdot \text{l-NME} \tag{1}$$

where:
- $W_0$ — endpoint energy in MeV
- $m_n - m_p \approx 1.293$ MeV — neutron-proton mass difference
- $E_C$ — Coulomb displacement energy (depends on nuclear charge distribution)
- l-NME — large vector NME (computable in NSM)

### Why is sNME important?

Despite being "small" (typically $< 0.1$ in absolute units), the sNME can influence β-spectral shapes and half-lives **quite strongly** because:

1. **Quadratic dependence:** The half-life $t_{1/2} \propto 1/(\text{l-NME} + \text{sNME})^2$, so even a small sNME shifts the total matrix element significantly.
2. **Small Q-value amplification:** For $^{99}\text{Tc}$, $[W_0 - (m_n - m_p)] = -0.487$ MeV — the negative mass difference nearly cancels the endpoint energy, making the sNME contribution relatively large.
3. **Shape dependence:** Different sNME values produce measurably different spectral shapes, especially below 100 keV.

## The CVC Value — An Idealization

The CVC prediction (Eq. 1) is an **idealization** — it pertains to a perfect many-body theory. In practice:
- NSM calculations with finite valence spaces give sNME $\approx 0$ (the valence space truncation kills it)
- The CVC value serves as a **reference**, not a physical prediction

## ESSM: sNME as a Fitting Parameter

### The Enhanced Spectrum Shape Method (ESSM)

Kumar et al. (2021) proposed using the sNME as an additional fitting parameter:

1. **Fix** $g_A$ at a trial value
2. **Vary** the sNME to match the **experimental partial half-life**
3. Because $t_{1/2} \propto (\text{l-NME} + \text{sNME})^{-2}$, the equation $t_{1/2}^{\text{theo}}(g_A, \text{sNME}) = t_{1/2}^{\text{exp}}$ typically yields **two solutions** for the sNME:
   - **sNME(c):** closer to the CVC value
   - **sNME(f):** further from the CVC value

### Two-Solution Structure

For each $g_A$ value, there are two sNME solutions (or none, if the experimental half-life cannot be reproduced). This creates two distinct spectral shape predictions:

- The **sNME(c)** solution is theoretically preferred (closer to CVC)
- The **sNME(f)** solution may also be physical — only experiment can decide

### Application to $^{99}\text{Tc}$

Song et al. (2025) applied ESSM to $^{99}\text{Tc}$ with three nuclear models:

| Model | $g_A^{\text{eff}}$ | sNME (c) [fm³] | sNME (f) [fm³] | $\chi^2_\nu$ | p-value |
|---|---|---|---|---|---|
| **CloseCVC (CVC fixed)** | 1.1 | 0.0674 | — | 1.106 | 0.12 |
| | 1.0 | 0.0681 | — | 1.244 | 0.004 |
| | 1.0 | 0.0651 | — | 1.644 | <0.0001 |
| **FarCVC (sNME fitted)** | 1.2 | −0.0698 | — | 1.066 | 0.22 |
| | 1.0 | −0.0669 | — | 1.097 | 0.14 |
| | 1.0 | −0.0683 | — | 1.095 | 0.14 |

Key findings:
- **Paulsen et al. (CloseCVC):** $g_A^{\text{eff}} = 1.526(92)$ — enhanced value
- **Song et al. (FarCVC/ESSM):** $g_A^{\text{eff}} \approx 1.0$–1.2 — quenched value
- The discrepancy arises from the **sNME treatment**, not from the data

## sNME in Other Transitions

Ramalho & Suhonen (2024) studied five second-forbidden nonunique decays:

| Transition | $g_A$ Sensitivity | sNME Sensitivity |
|---|---|---|
| $^{60}\text{Fe} \to {}^{60}\text{Co}$ | Strong | Weak |
| $^{94}\text{Nb} \to {}^{94}\text{Mo}$ | Strong (jj45pnb) | Weak |
| $^{98}\text{Tc} \to {}^{98}\text{Ru}$ | Strong | Moderate (below 100 keV) |
| $^{126}\text{Sn} \to {}^{126}\text{Sb}$ | Weak | **Strong** |
| $^{129}\text{I} \to {}^{129}\text{Xe}$ | Weak | **Strong** |

For $^{99}\text{Tc}$, the $g_A$ sensitivity is **strong** and the sNME sensitivity is **moderate** — making it an excellent candidate for combined $g_A$ + sNME fitting.

## References

| Ref | Citation | Use |
|---|---|---|
| [ESSM] | Kumar, Srivastava, Suhonen, *Eur. Phys. J. A* **57**, 225 (2021) | ESSM formalism |
| [Song2025] | Song et al., arXiv:2506.17544 (2025) | $^{99}$Tc with ESSM |
| [Ramalho2024] | Ramalho, Suhonen, *Front. Phys.* **12**, 1455778 (2024) | sNME in five second-forbidden decays |
| [Behrens] | Behrens & Bühring, *Electron Radial Wave Functions and Nuclear Beta Decay* (1982) | sNME formalism |
