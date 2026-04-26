---
title: "03 – Finite Size Effects"
date: 2026-04-26
tags:
  - beta-spectrum/correction/finite-size
status: active
aliases: [Finite Nuclear Size, L_0, U, DFS]
cssclasses: []
related_components: [[01-fermi-function]], [[10-nuclear-structure]]
---

# Finite Nuclear Size Effects — $L_0$, $U(Z,W)$, and $DFS$

## Overview

The Fermi function ([[01-fermi-function]]) assumes a point charge nucleus. Real nuclei have finite extent and possibly deformation — these corrections account for the difference in electric potential inside the nuclear volume.

> [!note] Two types of "finite size"
> 1. **Electrostatic** ($L_0$, $U$, $DFS$): Potential difference when moving from point charge to realistic nuclear shape — applied *after* evaluating wave function at origin
> 2. **Convolution**: Integration of leptonic wave functions over the nuclear volume — appears in shape factor $C(Z,W)$ ([[10-nuclear-structure]])

## Nuclear Radius

$$R = \sqrt{5/3} \, \langle r^2 \rangle^{1/2} \tag{14}$$

where $\langle r^2 \rangle^{1/2}$ is the root-mean-square charge radius from electron scattering experiments.

---

## $L_0(Z,W)$ — Uniformly Charged Sphere

### Low-Z Expansion (Behrens & Bühring)

$$L_0 \simeq 1 \mp \alpha Z W R + (\alpha Z)^2 - \frac{13}{60} - \frac{\alpha Z R}{2W} \tag{15}$$

### Wilkinson's Analytical Parametrization (Eq. 16)

Accurate to $10^{-4}$ for $p \leq 45$, $|Z| \leq 60$:

$$L_0(Z,W) = 1 + \frac{13}{60}(\alpha Z)^2 \mp \frac{\alpha Z W R(41 - 26\gamma)}{[15(2\gamma-1)]} $$
$$\mp \frac{\alpha Z R \gamma(17 - 2\gamma)}{[30W(2\gamma-1)]} + \frac{R}{W}\sum_{n=0}^5 a_n (WR)^n + 0.41(R-0.0164)(\alpha Z)^{4.5}$$

with $a_n = \sum_{x=1}^6 b_{x,n} (\alpha Z)^x$ (Eq. 17). Coefficients in **Tables I** (electrons) and **II** (positrons); signs of odd powers flipped for positrons.

> [!warning] Magnitude
> The effect ranges from a few 0.1% up to several percent — everywhere highly significant. Cannot be neglected for precision work at $10^{-4}$ level.

### Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/finite_size.py` → `FiniteSizeL0(Z, A)` class
> 
> Uses Wilkinson parametrization (Eq. 16) with coefficients from Hayen Tables I-II. The nuclear radius is computed via `utils.nuclear_radius(A)`.

---

## $U(Z,W)$ — More Realistic Charge Distributions

The uniform sphere approximation is insufficient for high precision. Two common models:

### Modified Gaussian (MG)
$$\rho_{MG}(r) = N_0\left[1 + A\left(\frac{r}{a}\right)^2\right]e^{-(r/a)^2} \tag{18}$$

with $N_0$ and $a$ determined by normalization and RMS radius (Eqs. 19–20).

### Fermi Distribution
$$\rho_F(r) \propto [1 + \exp((r-c)/a)]^{-1}, \quad a \simeq 0.55\text{ fm}$$

The correction is defined as ratio:

$$\frac{U(Z,W)}{L_0} = \frac{\alpha_{02}' + \alpha_{-1}'^2}{\alpha_{1}^2 + \alpha_{-1}^2} \tag{23}$$

### Analytical Approximation for Low Z (Eq. 25)

$$U(Z,W) \approx 1 + \alpha Z W R \Delta_1 + \frac{\gamma}{W}\alpha Z R \Delta_2 + (\alpha Z)^2 \Delta_3 - (WR)^2 \Delta_4$$

### Wilkinson's Fermi Distribution Fit (Eqs. 29–30)

$$U(Z,W) = 1 + \sum_{n=0}^\infty a_n p^n$$

with $a_0, a_1, a_2$ given as quadratic polynomials in $Z$.

> [!tip] Magnitude
> The effect of $U(Z,W)$ is ~0.1% for medium-high masses — cannot be neglected at the $10^{-4}$ precision level.

### Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/finite_size.py` → `ChargeDistributionU(Z, A)` class
> 
> Implements Eq. (25) with $\Delta_i$ coefficients computed from nuclear model parameters.

---

## $DFS(Z,W,\beta^2)$ — Nuclear Deformation

For axially symmetric deformations:

$$R(\theta,\phi) = R_0\left[1 + \beta_2 Y_{20}(\theta,\phi)\right] \tag{32}$$

with $\beta_2 > 0$ prolate, $\beta_2 < 0$ oblate. The intrinsic quadrupole moment:

$$Q_0 = \frac{\alpha}{\pi}\sqrt{\frac{5}{3}} R Z \beta_2(1 + 0.16\beta_2) \tag{33}$$

The correction is a ratio of angle-averaged potentials (Eq. 40):

$$DFS(Z,W,\beta) = \frac{L_0^*(Z,W)}{L_0(Z,W)} \tag{40}$$

where $L_0^*$ integrates over a continuous superposition of uniformly charged spheres (Eq. 39).

> [!tip] Magnitude
> Effect can reach several parts in $10^3$. Partially canceled by nuclear structure deformation correction $D_C$ ([[10-nuclear-structure]]). Figure 1 of Hayen shows energy dependence reversal between β⁻ and β⁺ for $\beta_2 = 0.2$.

### Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/finite_size.py` → includes DFS in the finite size module (combined with L₀)
> 
> The deformation parameter $\beta_2$ is typically provided by nuclear data tables or shell model calculations.

---

## Summary of Corrections Flow

```
Fermi function F₀(Z,W)   [point charge]
        ↓
    × L₀(Z,W)            [uniform sphere correction]
        ↓
    × U(Z,W)             [realistic charge distribution]
        ↓
    × DFS(Z,W,β²)         [nuclear deformation]
```

## References in Hayen Paper

- Section IV.A.1 (Eq. 15–17): $L_0$ uniform sphere
- Section IV.A.2 (Eq. 18–30): $U$ charge distributions  
- Section IV.A.3 (Eq. 32–40): $DFS$ deformation
