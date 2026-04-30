---
title: "06 – Radiative Corrections"
date: 2026-04-26
tags:
  - beta-spectrum/correction/radiative
status: active
aliases: [Outer Radiative, δ_R, Sirlin Function]
cssclasses: []
related_components: [[10-nuclear-structure]], [[07-recoil-effects]]
---

# Outer Radiative Corrections — $R(W,W_0) = 1 + \delta_R(W,W_0)$

## Overview: Inner vs. Outer Corrections

The radiatively corrected spectrum (Hayen et al., Eq. 46):

$$\frac{d\Gamma}{dW} = \frac{d\Gamma_0}{dW}(1 + \Delta_R^{V/A})(1 + \delta_R(W,W_0))$$

| Type | Symbol | Dependence | Treatment |
|---|---|---|---|
| **Inner** | $\Delta_R$ | Nuclear-independent, absorbed into effective coupling constants | $\Delta_V^R = (2.361 \pm 0.038)\%$ (Marciano & Sirlin 2006) — NOT part of spectrum shape calculation |
| **Outer** | $\delta_R(W,W_0)$ | Energy-dependent, nucleus-dependent | This is the correction implemented in code ([[06-radiative-corrections]]) |

## Total Spectral Influence (Eq. 47–48)

$$R(W,W_0) = 1 + \frac{\alpha}{2\pi}\left[g(W_0,W) - 3\ln\left(\frac{m_p}{2W_0}\right)\right] + Z\alpha^2\sum_i \Delta_i(W) + Z^2\alpha^3 \delta_{3h} + ... \tag{48}$$

### Order α — Sirlin's $g(W_0, W)$ Function (Eq. 50–52)

The universal correction (same for electrons and positrons, Fermi and GT):

$$g(W_0,W) = 3\ln(m_p) - \frac{3}{4} + \frac{4\beta}{1+\beta}\left[4\tanh^{-1}\beta - 1\right] $$
$$+ \frac{W_0-W}{\beta}\left(-\frac{3}{3W} + \ln[2(W_0-W)]\right) $$
$$- \frac{(W_0-W)^2}{6W^2} + 4\frac{\tanh^{-1}\beta}{\beta}(1+\beta^2) - 4\tanh^{-1}\beta$$

where $\beta = p/W$.

The Spence function (dilogarithm):
$$L_s(x) = \int_0^x \frac{\ln(1-t)}{t} dt = -Li_2(x) = -\sum_{k=1}^\infty \frac{x^k}{k^2} \tag{51}$$

Large $W_0$ limit:
$$g(W_0 \to \infty, W) = 3\ln\left(\frac{M_N}{2W_0}\right) + \frac{81}{10} - \frac{4\pi^2}{3} \tag{52}$$

> [!warning] Endpoint divergence
> $g(W_0, W)$ has a **logarithmic divergence at $W = W_0$** related to soft real photon emission. This must be handled via resummation.

### Soft Photon Resummation (Eq. 53–54)

To handle the endpoint divergence:
$$t(\beta)\ln(W_0 - W) \to (W_0-W)^{t(\beta)} - 1 \tag{53}$$

where:
$$t(\beta) = \frac{2\alpha}{\pi}\left[\frac{\tanh^{-1}\beta}{\beta} - 1\right] \tag{54}$$

For tritium the effect is negligible (low endpoint). For higher energies, the $ft$ value shift can reach several parts in $10^4$.

### Order Zα² Correction (Eq. 55–62)

Three dominant Feynman diagrams: vacuum polarization of bremsstrahlung photon and two electron propagator renormalizations. The axial vector component contributes specifically to $\Delta_3$:

$$\delta_2(Z,W) = Z\alpha^2 \sum_{i=1}^4 \Delta_i(W) \tag{55}$$

The leading term $\Delta_1$ depends on nuclear form factor $F(q^2)$ and is split into:
- **Nucleus-independent** part $\Delta_1^0$: energy-dependent terms combined with $\Delta_4$
- **Nuclear-dependent** part $\Delta_1^F$: explicit integral forms for uniformly charged sphere

For modified Gaussian model, differences between models for superallowed Fermi transitions ($^{10}$O to $^{54}$Co) are < $10^{-4}$.

### Order Z²α³ — Sirlin's Heuristic Formula (Eq. 63–67)

$$\delta_{3h} = Z^2\alpha^3 \left[a\ln\left(\frac{\Lambda}{W}\right) + bf(W) + \frac{4\pi}{3}g(W) - 0.649\ln(2W_0)\right]$$

with $a = \frac{\pi}{3} - \frac{\sqrt{3}}{2}$, $b = \frac{4}{3\pi^4}(\pi^2 - \gamma_E) - \frac{\pi^2}{18}$.

### Higher Order (Eq. 68–69)

$$\delta_{Z^n\alpha^{n+1}} \approx Z^n\alpha^{n+1} K_{nm} \ln^{m-n}\left(\frac{\Lambda}{W}\right), \quad K_{nm} \approx 0.5 \tag{68}$$
$$\sum_{n=3}^\infty Z^n\alpha^{n+1} = \frac{Z^3\alpha^4}{1-Z\alpha} \tag{69}$$

Relevant for higher $Z$ nuclei only.

## Neutrino Radiative Corrections (Eq. 70–71)

Even though neutrinos have no direct electromagnetic interaction, they receive indirect corrections via virtual photon exchange and energy conservation from inner bremsstrahlung:

$$R_\nu = 1 + \frac{\alpha}{2\pi}\left[3\ln\left(\frac{m_p}{m_e}\right) + \frac{23}{4} - Li_2\left(-\frac{8}{\hat{\beta}}\right) + ...\right]$$

where $\hat{W} = W_0 - W_\nu$. Much smaller than electron spectrum correction — relevant for reactor neutrino oscillation studies.

## Internal Bremsstrahlung (Eq. 72–74)

An additional photon in the final state: $n \to p + e^- + \bar{\nu}_e + \gamma$

Photon ejection probability per unit energy:
$$\Phi(W,\omega) = \frac{\alpha p}{\pi\omega p'}\left[\frac{W+W'^2}{p^2}\ln\left(\frac{W+p}{W'-p}\right) - 2\right] \tag{73}$$

Experimental branching ratios: $BR_{\beta\gamma} \sim 3\times 10^{-3}$ (neutron), $\sim 2\times 10^{-3}$ ($^{32}$P).

## Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/radiative.py` → `RadiativeCorrection(W0, Z, A, use_endpoint_resummation)` class
> 
> Implements outer radiative corrections per Sirlin (1967) and Sirlin (1987):
> - **O(α)**: Universal Sirlin function $g(W,W_0)$ with `scipy.special.spence` for dilogarithm
> - **O(Zα²)**: Finite nuclear size correction $\delta_2(Z,W)$ using modified Gaussian model
> - Soft-photon resummation at endpoint via Sirlin (1987) prescription
> - Small-beta Taylor expansion for numerical stability at threshold
> - Optional endpoint resummation flag

### API Usage

```python
from beta_spectrum import RadiativeCorrection

# Light nucleus (tritium-like)
rad = RadiativeCorrection(
    W0=1.02,
    Z=1,
    use_endpoint_resummation=True
)

# Heavy nucleus (uranium) with nuclear model correction
rad = RadiativeCorrection(
    W0=5.0,
    Z=92,
    A=238,  # affects Delta_F nuclear model term
    use_endpoint_resummation=True
)

values = rad(W_grid)
```

## Magnitude Summary

| Term | Typical Size | Notes |
|---|---|---|
| $\alpha/2\pi \cdot g$ | ~1–3% | Universal, dominant contribution |
| $Z\alpha^2 \Delta_i$ | ~0.01–0.1% | Nuclear-form-factor dependent |
| $Z^2\alpha^3 \delta_{3h}$ | < 0.01% | Only for heavy nuclei |
