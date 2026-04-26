---
title: Beta Spectrum Wiki
date: 2026-04-26
tags:
  - vault/overview
status: active
---

# Beta Spectrum Calculator — Knowledge Base

This Obsidian vault serves as the scientific reference and development workspace for the **beta-spectrum calculator** project (`beta_shape_pnpi`).

## Purpose

Organize the physics theory behind allowed beta decay spectrum calculations into structured, interconnected notes that directly inform the implementation of correction components in the codebase.

## Source Material

| Reference | Type | Location |
|---|---|---|
| Hayen et al., RevModPhys **90**, 015008 (2017) | Review article | `pdf_source/Hayen2017.pdf` / [[Hayen2017_summary]] |

## Notes Structure

All notes are flat in this directory, linked via wikilinks.

### Core Theory

- [[00-beta-spectrum-overview]] — Master equation and overview of all corrections
- [[02-phase-space]] — Baseline spectral shape $pW(W_0-W)^2$

### Correction Components (mapped to calculator modules)

| # | Note | Calculator Module | Status |
|---|---|---|---|
| 1 | [[01-fermi-function]] | `components/fermi.py` | ✓ implemented |
| 2 | [[03-finite-size]] | `components/finite_size.py` | ✓ implemented |
| 3 | [[04-screening-correction]] | `components/screening.py` | ✓ implemented |
| 4 | [[05-exchange-correction]] | `components/exchange.py` | ✓ implemented |
| 5 | [[06-radiative-corrections]] | `components/radiative.py` | ✓ implemented |
| 6 | [[07-recoil-effects]] | `components/recoil.py` | ✗ not yet implemented |

### Nuclear Structure & Atomic Effects

- [[10-nuclear-structure]] — Shape factors, form factors, impulse approximation
- [[11-isospin-breakdown]] — Shake-off/shake-up processes, endpoint shift
- [[12-atomic-overlap]] — Bahcall correction, bound-state β decay
- [[13-chemical-effects]] — Molecular environment influence on spectrum

### References

- [[refs/hayen2017]] — Metadata and key equations from Hayen et al. (2017)

## Calculator Architecture

```
SpectrumComponent (ABC, base.py)
 ├── PhaseSpace              → phase_space.py
 ├── FermiFunction           → fermi.py
 ├── FiniteSizeL0            → finite_size.py
 ├── ChargeDistributionU     → finite_size.py
 ├── ScreeningCorrection     → screening.py
 ├── ExchangeCorrection      → exchange.py
 └── RadiativeCorrection     → radiative.py

BetaSpectrum (spectrum.py) — multiplicative composition + from_config()
BetaSpectrumAnalyzer (spectrum.py) — plotting, CSV export
```

## Key Variables and Units

- $W$ — total electron energy in units of $m_e c^2$ (natural units, $m_e = 1$)
- $T$ or $E_{\text{kin}}$ — kinetic energy (keV or MeV)
- $W_0$ — endpoint total energy
- $Z_{\text{parent}}$, $Z_{\text{daughter}}$ — atomic numbers
- $\alpha \approx 1/137.036$ — fine structure constant
- $\gamma = \sqrt{1 - (\alpha Z)^2}$

Conversion: $W = 1 + T / m_e c^2$, with $m_e c^2 = 510.998950$ keV.
