---
title: "01 – Fermi Function"
date: 2026-04-26
tags:
  - beta-spectrum/correction/fermi
status: active
aliases: [Coulomb Correction, F_0]
cssclasses: []
related_components: [[03-finite-size]], [[06-radiative-corrections]]
---

# Fermi Function — Coulomb Interaction $F_0(Z,W)$

## The Point Charge Fermi Function

The transition matrix element with all-order Coulomb interactions resummed (Hayen et al., Eq. 5):

$$\mathcal{M}_{fi} = -2\pi i \delta(E_f - E_i) \langle f | T e^{-i\int_0^\infty dt H_Z(t)} H_\beta(0) T e^{-i\int_{-\infty}^0 dt H_{Z'}(t)} | i \rangle$$

Keeping only the final-state Coulomb interaction (daughter nucleus charge $Z$):

$$F_0(Z,W) = 4(2pR)^{2(\gamma-1)} e^{\pi y} \frac{|\Gamma(\gamma + iy)|^2}{(\Gamma(1+2\gamma))^2} \tag{6}$$

with:
$$\gamma = \sqrt{1 - (\alpha Z)^2}, \quad y = \pm \alpha Z W / p$$

Here $R$ is the cutoff radius representing daughter nucleus radius — necessitated by divergence at origin for point charge. Most of its dependence is absorbed into corrections ([[03-finite-size]]).

## Coulomb Amplitudes

Expanding radial Dirac wave functions near the origin (Eq. 12):

$$\left\{\frac{f_\kappa(r)}{g_\kappa(r)}\right\} = \alpha_\kappa \{(2|\kappa|-1)!!\}^{-1}(pr)^{|\kappa|-1} \sum_{n=0}^\infty \frac{a_{\kappa n}}{b_{\kappa n}} r^n$$

The Coulomb information is encoded in **Coulomb amplitudes** $\alpha_\kappa$. The dominant electron component at the origin:

$$F_0 L_0 = \frac{\alpha_{-1} + \alpha_1^2}{2p^2} \tag{13}$$

This is valid regardless of nuclear charge distribution shape.

## Alternative Definitions

| Convention | Prefactor for $L_0$ | Style |
|---|---|---|
| **Wilkinson / Blatt & Weisskopf** | $(1+\gamma)/2$ | Absorbs into $L_0$ |
| **Behrens-Bühring / Schopper** (used in calculator) | 4 | All corrections separate |

The error from neglecting the $(1+\gamma)/2$ prefactor is ~0.5% for $Z=20$.

## Numerical Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/fermi.py` → `FermiFunction(Z, A)` class
> 
> Uses `scipy.special.loggamma` for numerical stability (avoids overflow of $\Gamma$ function). The implementation follows Eq. (6) with log-space evaluation:
> ```python
> log_F0 = ln(4) + 2*(γ-1)*ln(2*p*R) + π*y 
>          + 2*Re[loggamma(γ+iy)] - 2*loggamma(1+2γ)
> F0 = exp(log_F0)
> ```

## Energy Dependence

$F_0$ is a strong function of $Z$:
- **β⁻ decay** ($y < 0$): wave functions attracted to nucleus → enhancement at low energy
- **β⁺ decay** ($y > 0$): repulsion → suppression at low energy, grows with energy

The effect ranges from ~1% for light nuclei to tens of percent for heavy ones.

## Related Corrections

After the point charge Fermi function is computed:
- [[03-finite-size]] corrects for non-pointlike nuclear charge ($L_0$, $U$, $DFS$)
- [[06-radiative-corrections]] adds QED radiative effects ($R = 1 + \delta_R$)
