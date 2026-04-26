---
title: "10 – Nuclear Structure Effects"
date: 2026-04-26
tags:
  - beta-spectrum/nuclear-structure
status: active
aliases: [Shape Factor, C(Z,W), Form Factors]
cssclasses: []
related_components: [[07-recoil-effects]], [[11-isospin-breakdown]]
---

# Nuclear Structure Effects — The Shape Factor $C(Z,W)$

## Overview

Nuclear structure effects are encapsulated in the shape factor $C(Z,W)$ which multiplies the phase space × Fermi function. This is distinct from electrostatic finite size corrections ([[03-finite-size]]) — it describes how nuclear form factors and matrix elements modify the spectrum beyond simple point-nucleus Coulomb interaction.

$$N(W)\,dW \propto p W (W_0 - W)^2 F(Z,W) C(Z,W) \tag{79}$$

## Two Formalisms Compared

| Approach | Style | Strengths |
|---|---|---|
| **Behrens-Bühring** (BB) | Couples leptonic + nuclear parts; rigorous expansion in $r^2$, $(m_eR)^a$, $(WR)^b$, $(\alpha Z)^c$ with tabulated functions $I(k_e,m,n,\rho)$ | Most precise, established computational machinery |
| **Holstein** (HS) | Decouples leptonic + nuclear parts; manifestly covariant form factors | Clean symmetry properties, intuitive classification |

Hayen et al. adopt BB's computational approach but present results in Holstein notation for clarity. See Appendix E of Hayen paper for detailed comparison.

## Multipole Expansion (Eq. 78)

Form factors $F_{KLs}(q^2)$ describe nuclear structure model-independently:

$$\langle f | V_\mu - A_\mu | i \rangle \propto \sum_{KMs} \sum_{L=K-1}^{K+1} (-1)^{J_f-M_f+M} (-i)^L F_{KLs}(q^2) $$
$$\times \sqrt{4\pi(2J_i+1)} \begin{pmatrix} J_f & K & J_i \\ -M_f & M & M_i \end{pmatrix} \frac{(qR)^L}{(2L+1)!!}$$

Since $(qR)^2$ is very small for beta decay, expansion typically stops after first order.

## Impulse Approximation (Eq. 80–81)

Nucleons treated as non-interacting particles coupling to weak vertex as free particles:

$$O_{KLs} = \sum_{\alpha,\beta} \langle \alpha | O_{KLs} | \beta \rangle a^\dagger_\alpha a_\beta \tag{80}$$
$$\langle f | O_{KLs} | i \rangle = \sum_{\alpha,\beta} \langle \alpha | O_{KLs} | \beta \rangle \langle f | a^\dagger_\alpha a_\beta | i \rangle \tag{81}$$

The one-body density matrix elements $\langle f | a^\dagger_\alpha a_\beta | i \rangle$ are calculable via shell model methods.

> [!warning] Limitations
> The impulse approximation neglects many-body effects:
> - **Meson exchange**: long history (Blin-Stoyle, Chemtob & Rho)
> - **Core polarization**: nuclear medium modifications
> 
> For allowed decays, results largely hold due to cancellations between core-polarization and meson exchange (Morita 1985), but individual effects can reach ~40%.

## Induced Currents — Form Factor Classification

See [[07-recoil-effects]] for detailed equations. Key points:

### First-Class Currents (CVC/PCAC)
- **Vector** $g_V$: conserved, from EM data via CVC
- **Weak magnetism** $(g_M - g_V)/2M$: determined by CVC from M1 transitions
- **Axial-vector** $g_A$: measured in neutron decay
- **Induced pseudoscalar** $g_P/2M$: PCAC predicts ~−229, quenched ~80% in nuclei

### Second-Class Currents (G-parity forbidden)
- **Induced scalar** $g_S/2M$: excluded by CVC for Fermi transitions
- **Induced tensor** $g_{II}/2M$: experimental null results; reduces to first-class component only

## Nuclear Matrix Element Operators (Table V of Hayen)

| Operator | Symbol | Type | $\Delta J$ |
|---|---|---|---|
| $M_F$ | Fermi matrix element | Vector | 0 |
| $M_{GT}$ | Gamow-Teller matrix element | Axial-vector | 0,1 (not 0→0) |
| $M_L$ | Longitudinal vector | — | depends on K |
| $M_{T}, M_{TL}$ | Transverse multipoles | — | depends on K |

## Validity of Impulse Approximation

- **Vector current**: conserved via CVC → impulse approximation reliable
- **Axial current**: broken by PCAC → quenching effects significant (up to 80% for $g_P$)
- Mesonic degrees of freedom absorbed into nuclear potential (Wilkinson 1974)

## Nuclear Deformation Correction $D_C(Z,W,\beta^2)$

The nuclear structure analogue of the electrostatic deformation correction ([[03-finite-size]]):
- Accounts for non-spherical nuclear shape in convolution terms
- Partially cancels with $DFS$ correction (electrostatic)
- Relevant for deformed nuclei ($\beta_2 \neq 0$)

## Calculator Status

> [!info] Not yet implemented
> The shape factor $C(Z,W)$ is NOT part of the current calculator. This note documents the theory for future implementation. For superallowed Fermi transitions, $C(Z,W) = |M_F|^2 \approx 1$ (isospin symmetry), so it can be approximated as unity in many cases.
