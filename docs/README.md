---
title: Beta Spectrum Calculator — Knowledge Base
date: 2026-04-26
tags:
  - vault/overview
status: active
---

# Beta Spectrum Calculator — Knowledge Base

This Obsidian vault serves as the scientific reference and development workspace for the **beta-spectrum calculator** project (`beta_shape_pnpi`).

## Purpose

Organize the physics theory behind β-decay spectrum calculations into structured, interconnected notes that inform the implementation of correction components in the codebase.

## Source Material

| Reference | Type | Location |
|---|---|---|
| Hayen et al., *Rev. Mod. Phys.* **90**, 015008 (2017) | Review article | [[refs/electromagnetic/2017_Hayen_analytical_beta_shape/summary]] |

## Directory Structure

```
docs/
├── README.md                          ← You are here
├── 00-beta-spectrum-overview.md       ← Master equation & full correction list
├── 02-phase-space.md                  ← Baseline spectral shape
├── corrections/                       ← All correction components
│   ├── 01-fermi-function.md           ← Fermi function F₀(Z,W)
│   ├── 02-finite-size.md              ← Nuclear size: L₀, U, DFS
│   ├── 03-screening.md                ← Atomic screening S(Z,W)
│   ├── 04-exchange.md                 ← Atomic exchange X(Z,W)
│   ├── 05-radiative.md                ← Outer radiative R(W,W₀)
│   ├── 06-recoil.md                   ← Recoil: R_N, Q, weak magnetism
│   ├── 07-detector-response.md        ← Detector convolution
│   ├── 08-atomic-overlap.md           ← Atomic mismatch r(Z,W)
│   └── 09-chemical-effects.md         ← Molecular environment
├── nuclear-structure/                 ← Nuclear theory
│   ├── index.md                       ← Overview
│   └── 01-shape-factors.md            ← C(Z,W), NMEs, impulse approximation
├── gA-extraction/                     ← gA extraction (new, in development)
│   ├── 01-ssm-algorithm.md            ← SSM step-by-step algorithm
│   ├── 02-ssnme.md                    ← sNME and ESSM formalism
│   ├── 03-tc99-case-study.md          ← 99Tc literature results
│   └── 04-implementation-todo.md      ← Implementation TODO list
└── refs/                              ← Paper summaries
    ├── electromagnetic/               ← CVC, radiative corrections
    ├── forbidden-decay/               ← Forbidden decay studies
    └── (organized by topic)
```

## Calculator Architecture

```
SpectrumComponent (ABC, base.py)
 ├── PhaseSpace              → 02-phase-space.md
 ├── FermiFunction           → corrections/01-fermi-function.md
 ├── FiniteSizeL0            → corrections/02-finite-size.md
 ├── ChargeDistributionU     → corrections/02-finite-size.md
 ├── ScreeningCorrection     → corrections/03-screening.md
 ├── ExchangeCorrection      → corrections/04-exchange.md
 ├── RadiativeCorrection     → corrections/05-radiative.md
 └── (recoil.py)             → corrections/06-recoil.md (not yet implemented)

BetaSpectrum (spectrum.py) — multiplicative composition + from_config()
BetaSpectrumAnalyzer (spectrum.py) — plotting, CSV export
```

## Correction Implementation Status

| # | Correction | Module | Note | Status |
|---|---|---|---|---|
| 1 | Phase space | `phase_space.py` | Baseline | ✓ implemented |
| 2 | Fermi function | `fermi.py` | Coulomb interaction | ✓ implemented |
| 3 | Finite size | `finite_size.py` | L₀, U, DFS | ✓ implemented |
| 4 | Screening | `screening.py` | Atomic screening | ✓ implemented |
| 5 | Exchange | `exchange.py` | Atomic exchange | ✓ implemented |
| 6 | Radiative | `radiative.py` | Outer RC | ✓ implemented |
| 7 | Recoil | `recoil.py` | R_N, Q, weak magnetism | ✗ not yet implemented |
| 8 | Shape factor | `shape_factor_gA.py` | C(Z,W) for forbidden | — planned |
| 9 | Detector | `detector_response.py` | Convolution | — planned |

## gA Extraction (New)

For forbidden nonunique decays (e.g., $^{99}\text{Tc}$), the shape factor depends on $g_A$ in a measurable way. See:

- [[gA-extraction/01-ssm-algorithm]] — The SSM algorithm for extracting $g_A$ from spectral shapes
- [[gA-extraction/04-implementation-todo]] — Step-by-step implementation plan

## Key Variables and Units

- $W$ — total electron energy in units of $m_e c^2$ (natural units, $m_e = 1$)
- $T$ or $E_{\text{kin}}$ — kinetic energy (keV or MeV)
- $W_0$ — endpoint total energy
- $Z_{\text{parent}}$, $Z_{\text{daughter}}$ — atomic numbers
- $\alpha \approx 1/137.036$ — fine structure constant
- $\gamma = \sqrt{1 - (\alpha Z)^2}$

Conversion: $W = 1 + T / m_e c^2$, with $m_e c^2 = 510.998950$ keV.

## Precision Goal

The calculator aims for **accurate to a few parts in $10^{-4}$ down to 1 keV** for low-to-medium $Z$ nuclei, following the approach of Hayen et al. (2017).
