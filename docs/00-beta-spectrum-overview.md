---
title: "00 – Beta Spectrum Overview"
date: 2026-04-26
tags:
  - beta-spectrum/overview
status: active
aliases:
  - Master Equation
cssclasses: []
---

# Beta Decay Spectrum — Full Expression

## The Master Equation

The fully differential allowed β spectrum shape, including all corrections needed for $10^{-4}$ precision (Hayen et al., Eq. 4):

$$N(W)\,dW = \frac{G_V^2 V_{ud}^2}{2\pi^3} \; K(Z,W,W_0,M) \; A(Z,W) \; C'(Z,W) \; p W (W_0 - W)^2 dW$$

or expanded:

$$N(W)\,dW = \frac{G_V^2 V_{ud}^2}{2\pi^3} F_0 L_0 U DFS R R_N Q S X r C D_C \; p W (W_0 - W)^2 dW \tag{4}$$

## Correction Factors

| Symbol | Name | Module / Note | Status |
|---|---|---|---|
| $F_0(Z,W)$ | Point charge Fermi function | [[01-fermi-function]] | ✓ |
| $L_0(Z,W)$ | Electrostatic finite size (uniform sphere) | [[03-finite-size]] | ✓ |
| $U(Z,W)$ | Charge distribution correction | [[03-finite-size]] | ✓ |
| $DFS(Z,W,\beta^2)$ | Nuclear deformation correction | [[03-finite-size]] | ✓ |
| $R(W,W_0) = 1 + \delta_R$ | Outer radiative corrections | [[06-radiative-corrections]] | ✓ |
| $R_N(W,W_0,M)$ | Recoil kinematic correction | [[07-recoil-effects]] | ✗ |
| $Q(Z,W,M)$ | Recoil electromagnetic correction | [[07-recoil-effects]] | — |
| $S(Z,W)$ | Atomic electron screening | [[04-screening-correction]] | ✓ |
| $X(Z,W)$ | Atomic exchange effect | [[05-exchange-correction]] | ✓ |
| $r(Z,W)$ | Atomic mismatch / overlap | [[12-atomic-overlap]] | — |
| $A(Z,W) = S \cdot X \cdot r$ | Combined atomic factors | [[12-atomic-overlap]] | — |
| $C(Z,W)$ | Nuclear structure shape factor | [[10-nuclear-structure]] | — |
| $D_C(Z,W,\beta^2)$ | Nuclear deformation (convolution) | [[10-nuclear-structure]] | — |

## Baseline: Phase Space

The uncorrected spectral shape is pure phase space ([[02-phase-space]]):

$$\frac{dN}{dW}\bigg|_{\text{PS}} = p W (W_0 - W)^2$$

with $p = \sqrt{W^2-1}$, $W$ total energy in units of $m_e$, and $W_0$ the endpoint.

## Precision Goal

The paper achieves **accurate to a few parts in $10^{-4}$ down to 1 keV** for low-to-medium $Z$ nuclei — extending previous work by nearly an order of magnitude [[Hayen2017_summary]].

> [!key] Why this matters
> A precision of ~0.5% is required when determining the Fierz term $b_{\text{Fierz}}$ to improve limits on BSM scalar/tensor couplings. Weak magnetism $b_{wm}$ must be included as experimental precision increases, otherwise it limits BSM sensitivity.

## Three Formalisms Compared

| Approach | Style | Strengths |
|---|---|---|
| **Behrens-Bühring** (BB) | Rigorous, numerical Dirac solutions | Most precise results |
| **Holstein** (HS) | Analytical, manifestly covariant | Clean symmetry properties |
| **Wilkinson** | Analytical parametrizations | Practical fits for dominant effects |

The Hayen paper combines BB's accuracy with HS's transparency.

## Conventions

- Natural units: $c = \hbar = m_e = 1$ throughout
- $Z$ is always positive; $\beta^-$ and $\beta^+$ distinguished via upper/lower signs ($\mp$)
- $W = E/m_ec^2 + 1$, $p = \sqrt{W^2 - 1}$
- $G_V$: vector coupling strength; $V_{ud} = \cos^2\theta_C$

> [!note] Two types of "finite size"
> The authors distinguish two corrections conflated under this term:
> 1. **Electrostatic** ($L_0, U, DFS$): Potential difference from point charge → realistic nuclear shape
> 2. **Convolution**: Integration of leptonic wave functions over the nuclear volume (appears in $C(Z,W)$)

## Energy Conversion

$$W = 1 + \frac{T}{m_ec^2}, \quad m_ec^2 = 510.998950\text{ keV}$$

In code: see `beta_spectrum/utils.py` → `T_to_W()`, `W_to_T()` functions.

## Calculator Implementation

Each correction factor $C_i(W)$ is a callable class inheriting from `SpectrumComponent`. The total spectrum is the product of all enabled components:

```
N(W) = PhaseSpace(W0) × FermiFunction(Z) × FiniteSizeL0(Z,A) × ... 
```

Configuration toggles are defined in `SpectrumConfig` (spectrum.py). See [[README]] for architecture details.
