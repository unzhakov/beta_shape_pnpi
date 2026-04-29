---
title: "02 – Phase Space"
date: 2026-04-26
tags:
  - beta-spectrum/correction/phase-space
status: active
aliases: [Phase Space, Baseline Spectrum]
cssclasses: []
related_components: [[01-fermi-function]], [[07-recoil-effects]]
---

# Phase Space — Baseline Spectral Shape $p W (W_0 - W)^2$

## The Unmodified Spectrum

The baseline β spectrum shape without any corrections is pure phase space:

$$\frac{dN}{dW}\bigg|_{\text{PS}} = p W (W_0 - W)^2 \tag{4, first factor}$$

where:
- $p = \sqrt{W^2 - 1}$ — electron momentum ($m_e c$ units)
- $W$ — total electron energy ($m_e c^2$ units)
- $W_0$ — endpoint (Q-value + $m_ec^2$)

## Neutrino Mass Effects

If the neutrino has non-zero mass $m_\nu$, the spectrum is modified near the endpoint:

$$p \to p_e = \sqrt{W^2 - 1}, \quad p_\nu = \sqrt{(W_0 - W)^2 - m_\nu^2}$$

The phase space factor becomes:
$$(W_0 - W)\sqrt{(W_0 - W)^2 - m_\nu^2}$$

This produces a characteristic distortion near $W = W_0 - m_\nu$, used in neutrino mass experiments (e.g., KATRIN with tritium β decay).

## Kurie Plot

The phase space subtracted spectrum gives the **Kurie function**:

$$K(W) = \sqrt{\frac{N(W)}{p W F(Z,W)}} \approx W_0 - W \quad (\text{linear for allowed transitions})$$

Deviations from linearity near the endpoint indicate:
- Non-zero neutrino mass
- Weak magnetism / recoil effects ([[07-recoil-effects]])
- New physics (Fierz term [[Hayen2017_summary#The Fierz Interference Term]])

## Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/phase_space.py` → `PhaseSpace(W0, transition_type)` class
> 
> Evaluates $p \times W \times (W_0 - W)^2$ multiplied by a **transition-type-dependent forbidden factor**. Supports optional neutrino mass.

## Transition Type and Forbidden Factors

The `PhaseSpace` calculator accepts a `transition_type` parameter that determines the multiplicative forbidden factor applied to the baseline phase space:

| `transition_type` | Forbidden Order | Forbidden Factor |
|---|---|---|
| `A`, `F1` | 0, 1 | $1$ |
| `F1U`, `F2` | 1, 2 | $p_\nu^2 + p_e^2$ |
| `F2U`, `F3` | 2, 3 | $p_\nu^4 + \frac{3}{10}p_\nu^2 p_e^2 + p_e^4$ |
| `F3U`, `F4` | 3, 4 | $p_\nu^6 + 7p_\nu^4 p_e^2 + 7p_\nu^2 p_e^4 + p_e^6$ |

The forbidden factor modifies the energy dependence of the spectrum, becoming significant at higher electron energies. See `[[10-nuclear-structure]]` for the full shape factor treatment.

## Energy Range

The phase space factor:
- Vanishes at $W = 1$ ($T = 0$, the threshold) — goes as $(W-1)^{1/2}$ via momentum
- Vanishes at $W = W_0$ (the endpoint) — quadratic behavior from neutrino energy
- Peaks roughly at $W \approx W_0 / 3$ to $W_0 / 2$ depending on $Z$

## Code API

```python
from beta_spectrum import PhaseSpace, T_to_W

ps = PhaseSpace(W0=T_to_W(5.0), transition_type="A")  # allowed transition, 5 MeV endpoint
values = ps(W_grid)                                    # returns p·W·(W₀-W)² × forbidden_factor
```
