---
title: "$^{99}$Tc Case Study — gA Extraction Results"
date: 2026-05-14
tags:
  - gA/extraction/tc99
  - beta-spectrum/forbidden/second
status: active
aliases: [Tc99 gA, 99Tc beta spectrum, SSM Tc99]
related_notes:
  - gA-extraction/01-ssm-algorithm
  - gA-extraction/02-ssnme
  - refs/2024_Paulsen_T99_beta_spectrum/summary
  - refs/2025_Song_tc99_measurement_for_g_v/summary
---

# $^{99}$Tc — The Premier Laboratory for gA Extraction

## Decay Properties

| Property | Value |
|---|---|
| Transition | $^{99}\text{Tc}(9/2_1^+) \to {}^{99}\text{Ru}(5/2_1^+)$ |
| Forbiddenness | Second-forbidden nonunique |
| Branching ratio | 99.99855(30)% |
| Q-value | 295.82(16) keV (Paulsen et al.) |
| Half-life | 211.5(1.1) $\times 10^3$ a (ENSDF) |
| Dominant transition | $2d_{5/2} \to 1g_{9/2}$ (OBTD $\approx 0.478$) |
| Sub-dominant | $1g_{7/2} \to 1g_{9/2}$ (OBTD $\approx 0.010$) |

## Two Contradictory Results

### Result 1: Paulsen et al. (2024) — Enhanced $g_A$

**Measurement:** MMC at PTB Berlin + LNHB Paris. Two independent setups, consistent results.

**Method:** Fixed sNME to CVC value (CloseCVC).

**Result:**
$$g_A^{\text{eff}} = 1.526(92)$$

This is an **enhanced** value — significantly above the free-nucleon value of 1.27.

**Follow-on values:**
- $g_V^{\text{eff}} = 0.376(5)$ — quenched (via half-life normalization)
- $g_A^{\text{eff, quenched}} = 0.574(36)$ — quenched
- $\bar{E}_\beta = 98.51(23)$ keV
- $\log ft = 12.3478(23)$

**Key claim:** This resolves an inconsistency — first-forbidden nonunique transitions show $g_A^{\text{eff}} \approx 0.9$, while fourth-forbidden show $g_A^{\text{eff}} \approx 0.9$. The $^{99}\text{Tc}$ result suggests $g_A$ is **enhanced** in second-forbidden decays.

### Result 2: Song et al. (2025) — Quenched $g_A$

**Measurement:** MMC at LLNL (KRISS sensor). Independent of Paulsen et al.

**Method:** ESSM — sNME fitted as free parameter (FarCVC).

**Result:**
$$g_A^{\text{eff}} \approx 1.0\text{--}1.2$$

This is a **quenched** value, consistent with first-forbidden and fourth-forbidden nonunique decays.

**Best fit:** glekpn Hamiltonian, $g_A^{\text{eff}} = 1.2$, sNME = −0.0698 fm³, $\chi^2_\nu = 1.066$, p-value = 0.22

**Key claim:** The Paulsen result is an artifact of fixing sNME to CVC. When sNME is fitted, the data prefer quenched $g_A$, consistent with the broader pattern.

## The Crucial Difference: sNME Treatment

| Aspect | Paulsen et al. | Song et al. |
|---|---|---|
| sNME treatment | Fixed to CVC value | Fitted as free parameter (ESSM) |
| $g_A^{\text{eff}}$ | 1.526(92) (enhanced) | 1.0–1.2 (quenched) |
| Nuclear models | NSM (GLEKPN) only | NSM (glekpn, jj45pnb) + MQPM |
| Half-life constraint | Separate normalization | Simultaneous shape + half-life fit |
| Consistency with 113Cd, 115In | Inconsistent | Consistent |

## What This Means for Your Experiment

The $^{99}\text{Tc}$ case demonstrates that **the treatment of sNME is the dominant theoretical systematic** in gA extraction. The two experiments (Paulsen and Song) measured essentially the same spectrum (Paulsen's data was used as input for Song's analysis) but obtained **contradictory results** because of different theoretical assumptions.

**For your data analysis:**
1. **Do not fix sNME to CVC** — treat it as a fitting parameter
2. **Use ESSM** — simultaneous fit to spectral shape + half-life
3. **Try multiple nuclear models** — glekpn, jj45pnb, and MQPM (if available)
4. **Report both sNME solutions** (c and f) for each model
5. **The "optimal" choice** (closest to CVC) is theoretically motivated but requires experimental validation

## Your Polynomial Fit: Connecting to SSM

You've fitted $C(W) = A \cdot W + B \cdot W^{-1} + C \cdot W^2$. This empirical parametrization can be **cross-checked** against the SSM prediction:

1. Compute $C_{\text{SSM}}(W; g_A)$ for a range of $g_A$ values
2. Fit each SSM curve to the polynomial form $aW + bW^{-1} + cW^2$ over your fit range
3. Extract the effective $(A, B, C)$ as a function of $g_A$
4. Compare your fitted $(A, B, C)$ to the $g_A$-dependent curve to identify which $g_A$ range is consistent

This provides a **model-independent bridge** between your empirical fit and the physical $g_A$ extraction.

## References

| Ref | Citation |
|---|---|
| [Paulsen2024] | Paulsen et al., *Phys. Rev. C* **110**, 055503 (2024) |
| [Song2025] | Song et al., arXiv:2506.17544 (2025) |
| [Kostensalo2017] | Kostensalo & Suhonen, *Phys. Rev. C* **96**, 024317 (2017) |
| [ESSM] | Kumar et al., *Eur. Phys. J. A* **57**, 225 (2021) |
