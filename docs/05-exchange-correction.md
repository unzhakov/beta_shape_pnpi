---
title: "05 – Exchange Correction"
date: 2026-04-26
tags:
  - beta-spectrum/correction/exchange
status: active
aliases: [Atomic Exchange, X(Z,W), Hayen-2018]
cssclasses: []
related_components: [[04-screening-correction]], [[12-atomic-overlap]]
---

# Atomic Electron Exchange Correction $X(Z,W)$

## Physical Origin

When a β electron is emitted from the nucleus, it can exchange with bound atomic electrons due to indistinguishability (Pauli exclusion principle). This **electron exchange** effect modifies the β spectrum at low energies where the continuum electron's wave function overlaps significantly with atomic orbitals.

> [!key] Magnitude
> The exchange contribution can exceed **20% in low-energy regions for higher Z** nuclei — it is one of the largest atomic corrections and absolutely critical for precision work below ~10 keV.

## Theory: Full Potential Approach

Exchange requires radial wave functions of both continuum and bound states integrated over **all space** (not just near origin). This is computationally demanding for analytical description because:
- Exchange needs integration over entire spatial domain
- Hydrogenic approach with simple screening is insufficient for high precision
- Must use realistic atomic potentials (Dirac-Hartree-Fock-Slater)

### Three Potential Models Compared (Hayen Fig. 9)

| Model | Description | Accuracy |
|---|---|---|
| **Unoptimized** | Basic Hartree-Fock-Slater, no exchange optimization | Lower bound |
| **Salvat** | Parametrized atomic potential | Intermediate |
| **Optimized** | Full potential with exchange optimization included | Target ($10^{-4}$) |

Screening strength is adjusted for best agreement with bound state energies (decreasing importance for higher $n$).

## Analytical Parametrisation — Hayen et al. 2018

The exchange correction $X(W)$ as a function of kinetic energy uses empirical fit coefficients from Table X of [[refs/hayen-2018]]. Coefficients are tabulated for all atoms with $Z = 2$ through $z = 120$.

### Shell Structure Effects

> [!tip] Figure 9 observation
> Clear shell closure effects observed in orbital exchange contributions: at noble gas Z values, the exchange correction shows characteristic kinks. These must be captured by interpolation between tabulated Z values.

## Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/exchange.py` → `ExchangeCorrection(Z)` class
> 
> Loads empirical fit coefficients from `data/exchange_coeff.csv` (119 rows, $Z = 2$–$120$). Uses interpolation between tabulated values.

### Data File

```
beta_spectrum/data/exchange_coeff.csv
├── Row count: 119 (Z=2 to Z=120)
├── Columns: Z, fit coefficients a₀...aₙ for energy parametrisation
└── Source: Hayen et al., Phys. Rev. D 97, 013006 (2018), Table X
```

### API Usage

```python
from beta_spectrum import ExchangeCorrection

exch = ExchangeCorrection(Z=90)  # Z of parent atom
values = exch(W_grid)           # returns correction factor at each W
```

## Energy Dependence

- **Dominant at very low energies** (T < ~10 keV): where continuum electron overlaps with atomic orbitals
- Decreases rapidly above ~50 keV (continuum wave function too short-ranged for significant overlap)
- For medium/high Z nuclei, shake-off probability in $s_{1/2}$ or $p_{1/2}$ state is < 0.1% — negligible

## Shake-Up Influence on Exchange

In low-Z nuclei (e.g., tritium):
- Shake-up probabilities can be significant (~25%)
- Exchange with higher-lying ns orbitals must be included even if unoccupied in initial state

For most practical cases, inclusion of shake-up changes the exchange correction at the few $10^{-4}$ level down to 0.5 keV for tritium only — otherwise negligible.

## Code Status: ✓ Implemented

The `ExchangeCorrection` class is fully implemented and tested. The CSV data file covers the full range from light (He, Z=2) to transuranic (U, Z=92+) elements.
