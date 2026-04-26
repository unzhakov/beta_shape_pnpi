---
title: "12 – Atomic Overlap Correction"
date: 2026-04-26
tags:
  - beta-spectrum/correction/atomic-overlap
status: active
aliases: [Atomic Mismatch, Bahcall Correction, r(Z,W)]
cssclasses: []
related_components: [[04-screening-correction]], [[13-chemical-effects]]
---

# Atomic Overlap / Mismatch Correction $r(Z,W)$

## Physical Origin

β decay results in a sudden change of nuclear potential ($Z \to Z\pm 1$). Initial and final atomic orbital wave functions only partially overlap — the atom must "rearrange" from its initial configuration to accommodate the new nuclear charge. This **atomic mismatch** (or shake-up) correction accounts for discrete excitations into higher bound states.

## Bahcall Correction (Eq. 166)

In first approximation, this reduces to a difference in atomic binding energies:

$$r(Z,W) = 1 - \frac{1}{W_0 - W} \frac{\partial^2 B(G)}{\partial Z^2} \tag{166}$$

where $B(G)$ is the total atomic binding energy for a neutral atom with $Z\pm 1$ protons. The second derivative relates to average excitation energy via:
$$\Delta E_{ex} = -\frac{1}{2}\frac{\partial^2 B(G)}{\partial Z^2}$$

### Parametrization (Eq. 169–170)

Using Desclaux (1973) numerical values:
$$K(Z) = -0.872 + 1.270 \, Z^{0.097} + 9.062\times 10^{-11} \, Z^{4.5} \tag{169}$$

The correction then becomes:
$$r(Z,W) = 1 - \frac{1}{W_0 - W^2}\left[B(G)'' + 2(C'_0 + C_1)\right] \tag{170}$$

## Magnitude Near Endpoint

> [!warning] Critical near endpoint for low-energy transitions
> Influence felt mainly near the endpoint where $\partial^2 B(G)/\partial Z^2$ can reach **a few hundreds of eV**.

| Isotope | Endpoint | Correction at 15 keV before endpoint |
|---|---:|---|
| $^{63}$Ni | 67.2 keV | ~**1%** |
| $^{241}$Pu | 20.8 keV | similarly important |
| Tritium ($^3$H) | 18.6 keV | relevant for KATRIN analysis |

The correction increases rapidly past the endpoint (where it's not physically meaningful — this signals breakdown of approximation).

## Bound State β Decay (Eq. 171)

Since the electron is created inside an electronic potential well, there exists a possibility for β decay to be captured into a bound atomic state:

$$\frac{\Gamma_b}{\Gamma_c} = \frac{\pi(\alpha Z)^3}{f(Z,W_0)} (W_0 - 1)^2 \Sigma \tag{171}$$

### Branching Ratios

| System | $\Gamma_b/\Gamma_c$ | Relevance for calculator |
|---|---:|---|
| Free neutron | $4.2\times 10^{-6}$ | Negligible |
| Tritium (T⁺) | ~**1%** | Relevant for precision neutrino mass experiments |
| Tritium (neutral T) | ~**0.5%** | — |
| Superallowed decays (MeV Q-value) | << $10^{-4}$ | Negligible |

> [!note] ft analysis only
> Bound state decay does **not** affect the β spectrum shape (it is a separate final state in S-matrix calculation). It enters when considering the **ft analysis**, so it's not directly relevant for `BetaSpectrum` but may matter for future calculator extensions.

## Calculator Implementation

> [!info] Not yet implemented
> The atomic mismatch correction $r(Z,W)$ is NOT currently part of the calculator. This note documents the theory for future implementation, particularly important for low-energy transitions (tritium, $^{63}$Ni, $^{241}$Pu).

## Related Notes

- [[04-screening-correction]] — Screening modification from shake-up processes
- [[11-isospin-breakdown]] — Shake-off process correction
- [[13-chemical-effects]] — Molecular environment influence on atomic overlap
