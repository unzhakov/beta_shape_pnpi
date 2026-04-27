---
title: "Hayen et al. (2017) — Reference Metadata"
date: 2026-04-26
tags:
  - beta-spectrum/references
status: active
aliases: [RRHP2017, RevModPhys.90.015008]
cssclasses: []
---

# Hayen et al. (2017) — "Beta-decay spectrum shape"

## Citation

> L. Hayen, K. Haseroth, N. Severijns, W. H. E. Schwarz, *Rev. Mod. Phys.* **90**, 015008 (2018).  
> DOI: [10.1103/RevModPhys.90.015008](https://doi.org/10.1103/RevModPhys.90.015008)

## PDF Location

`~/LBML_cloud/ОПЯД/beta_corrections/Refs/Hayen2017.pdf` (2.1 MB, 63 pages)

Summary: `docs/Hayen2017_summary.md`  
BibTeX references: `data/hayen2017_refs.bib` (319 entries extracted from PDF)

## Key Equations Reference

| Eq # | Section | Description | Wiki Note |
|------|---------|-------------|-----------|
| (1)–(5) | II | Dirac radial wave functions, Coulomb phase | [[01-fermi-function]] |
| (6) | II.2 | Fermi function $F_0(Z,W)$ | [[01-fermi-function]] |
| (7) | III | Phase space factor | [[02-phase-space]] |
| (9)–(15) | IV.A | Finite size: uniform sphere, Coulomb energy | [[03-finite-size]] |
| (16) | IV.A.1 | $L_0(Z,W)$ — finite size correction to Fermi function | [[03-finite-size]] |
| (17)–(24) | IV.A.2 | Modified Gaussian charge distribution, fit parameters | [[03-finite-size]] |
| (25), (29) | IV.B | $U(Z,W)$ — diffuse nuclear surface correction | [[03-finite-size]] |
| (40) | IV.C | $DFS$ — deformation correction to electrostatic potential | [[03-finite-size]] |
| (41)–(45) | V | Recoil kinematic corrections, weak magnetism terms | [[07-recoil-effects]] |
| (46)–(47) | VI.A | Radiative corrections $\delta_R'(W,W_0)$, Sirlin $g(W_0,W)$ | [[06-radiative-corrections]] |
| (53) | VI.A.2 | Soft-photon resummation $(W_0-W)^{t(\beta)} - 1$ | [[06-radiative-corrections]] |
| (78)–(81) | VII.A | Multipole expansion, form factors $F_{KLs}$, impulse approximation | [[10-nuclear-structure]] |
| (83), (92)–(104) | VI.E | Shape factor $C(Z,W)$ in Holstein notation | [[10-nuclear-structure]] |
| (107) | VI.G | $\Lambda$ factor for ft value analysis, weak magnetism contribution | [[07-recoil-effects]] |
| (113) | VI.F | Isovector correction $C_I(Z,W)$ — isospin symmetry breakdown | [[11-isospin-breakdown]] |
| (140), (144) | VII.A.2 | Atomic screening $S_{\text{Sal}}(Z,W)$ with Salvat potential | [[04-screening-correction]] |
| (157) | VII.E.3 | Analytical fit for exchange correction $X(Z,W)$ | [[05-exchange-correction]] |
| (160), (163), (164) | VII.C–D | Shake-up, shake-off processes | [[12-atomic-overlap]], [[11-isospin-breakdown]] |
| (166) | VII.D | Bahcall atomic mismatch correction $r(Z,W)$ | [[12-atomic-overlap]] |
| (171) | VII.E | Bound state β decay branching ratio $\Gamma_b/\Gamma_c$ | [[12-atomic-overlap]] |
| (176) | VII.F.3 | Molecular screening correction $\Delta S_{\text{Mol}}$ | [[13-chemical-effects]] |
| (178) | IX | Normalized shape factor $S(W)$ for BSM searches | [[07-recoil-effects]] |

## Table References

| Table | Content | Wiki Note |
|-------|---------|-----------|
| I | Leptonic radial wave function expansion coefficients | — |
| II | Leptonic and nuclear quantum numbers, multipole classification | [[10-nuclear-structure]] |
| III | Form factor mapping BB ↔ HS notation | [[10-nuclear-structure]] |
| IV | Numerical values for $L_0$ correction (Z=2–54) | [[03-finite-size]] |
| V | Nuclear matrix element operators ($M_F$, $M_{GT}$, etc.) | [[10-nuclear-structure]] |
| VI | Numerical values for screening correction $S(Z,W)$ at 67.2 keV (Ni) | [[04-screening-correction]] |
| VII | Complete overview of all β spectrum shape corrections | [[00-beta-spectrum-overview]] |

## Key Numbers to Remember

- **Target precision**: $10^{-4}$ down to 1 keV for low-to-medium Z nuclei
- **Exchange correction**: Can exceed **20%** at low energy for high Z (Table X, Eq. 157)
- **Radiative corrections**: Sirlin's $g(W_0,W)$ diverges logarithmically at endpoint → handled via soft-photon resummation
- **Induced pseudoscalar** $g_P(0)$: PCAC predicts ~−229 for free nucleon, quenched by up to 80% in nuclear medium
- **Weak magnetism uncertainty**: Can shift $f_A/f_V$ values by several parts in $10^3$

## Crosschecks in Paper

- **Superallowed Fermi decays** (§VIII): Excellent agreement with Towner & Hardy (2015) f values, all residuals in few $\times 10^{-4}$ region
- **Mirror decays**: Exquisite agreement for vector sector (< $4\times 10^{-4}$); axial vector depends on nuclear model choice

## Related Notes

- [[00-beta-spectrum-overview]] — Master equation and full corrections table (Table VII)
- [[10-nuclear-structure]] — Shape factor formalism, BB vs HS comparison (§E in paper)
- [[07-recoil-effects]] — Recoil kinematics, weak magnetism, induced currents

## Reference Extraction Tool

BibTeX was extracted using `extract_refs.py` (local regex mode), which parses pdftotext output to extract journal, volume, pages, arXiv IDs for 83% of references. For papers with non-standard citation formats, use `--llm` flag to invoke Qwen3.6 via llama-server API.
