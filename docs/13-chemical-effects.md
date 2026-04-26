---
title: "13 – Chemical Effects"
date: 2026-04-26
tags:
  - beta-spectrum/correction/molecular
status: active
aliases: [Molecular Screening, Molecular Exchange, Saenz-Froelich]
cssclasses: []
related_components: [[12-atomic-overlap]], [[04-screening-correction]]
---

# Chemical / Molecular Effects on β Spectrum

## Physical Origin

In many experiments (especially precision neutrino mass measurements), the decaying atom is bound within a molecule. This modifies:
- Electronic structure via molecular orbitals
- Coulombic final state interactions from additional electrons and spectator nuclei
- Energy transfer to molecular rovibrational states during decay

> [!warning] Systematic error impact
> In most precise analyses, this introduces a **systematic error of 0.05%**. Primarily relevant for tritium β-decay experiments (KATRIN).

## 1. Recoil Corrections in Molecules (§F.1)

When β decay occurs inside a molecule:
- The recoiling daughter nucleus moves through the molecular potential, not vacuum
- Potential described by Born-Oppenheimer energy curve, modelled with Lennard-Jones type potential
- Daughter nucleus kicks molecule into (predissociative) rovibrational state

**Dissociation channels:**
| Channel | Energy scale | Relevance |
|---|---:|---|
| Rovibrational states | ~few eV | Negligible at target precision |
| Electronic excitation → resonant continuum | tens–hundreds of eV | — |

Integrated probability depends on endpoint energy (higher endpoint → higher recoil). At $10^{-4}$ precision, this effect is neglected as a small correction on already small corrections.

For β-ν correlation measurements ($a_{\beta\nu}$):
- Dissociation probability decreases with increasing Z (recoil energy $\propto 1/m$)
- Can be partially included via Monte Carlo simulation in Lennard-Jones potential

## 2. Q Value Shift (§F.2)

Molecular effects add further decrease to Q value from **rotational and vibrational excitations**.

For $T_2$ (tritium molecule):
- Excitation possibilities change to a broad continuum with resonances
- Width of populated rovibrational spectrum: **few eV** — typically negligible
- Differences between tritium-substituted molecules: **few eV**, including atomic tritium

### Decay Rate Change (Eq. 172)

$$\frac{\Delta \lambda}{\lambda} \approx \frac{3 \Delta W_0}{W_0} \left( 1 + \gamma \frac{W_0}{6} \right) \tag{172}$$

where $\Delta W_0$ is difference in mean endpoint energy after averaging over all final states between two chemical states.

### Expected Error Level

Spectral shape depends on $W_0^2$, so relative error goes as $2 \times \sigma_Q / Q$. For ft values ($\propto Q^5$):
- Lowest energy transitions (tritium, $^{63}$Ni, $^{241}$Pu): **few $\times 10^{-4}$ to $10^{-3}$** level error for $\sigma_Q \sim$ tens of eV

## 3. Molecular Screening (§F.3) — Saenz & Froelich (1997a)

All Coulombic effects treated equally to first order, split into electronic and nuclear parts.

### Molecular Screening Correction (Eq. 176)

$$\Delta S_{\text{Mol}} = \alpha \left[ \frac{W}{p^2} + \frac{1}{W} \sum_S \frac{Z_S - Z_{S,\text{in}}}{R_{S,e}} - Z_{\text{eff}} \langle r^{-1} \rangle_{\text{Val}} \right]$$

where:
- $Z_S - Z_{S,\text{in}}$ = screened spectator nucleus charge
- $Z_{\text{eff}} = Z_{\text{Val}} - (Z - Z_{\text{in}})$
- $\langle r^{-1} \rangle_{\text{Val}}$ = average inverse distance of all valence electrons relative to decaying nucleus
- Sum over S runs through spectator nuclei, $R_{S,e}$ is equilibrium distance

### Example: $^{45}$Ca in CaCl₂ (Ca(II), doubly oxidized)
- Bond length: 2.437 Å → natural units $R_{\text{Cl},e} \sim 600$
- First term $\propto 3\times 10^{-5}$ — two orders of magnitude smaller than bare $V_0$
- Second term ($Z_{\text{eff}} = 2$): also ~**$2\times 10^{-5}$**, opposite sign
- Overall molecular deviation from atomic structure: **~$2\times 10^{-5}$**

Since screening is typically at (sub)percent level max, molecular deviations have upper limit of **few $\times 10^{-4}$**.

## 4. Molecular Exchange Effect (§F.4)

In a molecule, electronic phase space is greatly enlarged and perturbed relative to single atomic state. Electron density near decaying atom can both increase (bond region) or decrease, modifying the overlap integral.

### LCAO Approximation

Molecular orbitals constructed from linear combination of atomic orbitals. Key cases:
- **Ionic bond** (valence electron nearly fully removed): exchange effect ≈ 0 (internuclear distances >> β Compton wavelength where wave function oscillates rapidly)
- **Reverse extreme** (valence electron density doubled): naively doubles the exchange correction; conservative treatment uses **100% error bar on last orbital contribution**

### Error Estimate

| Isotope | Low-energy effect | High-energy behavior |
|---|---:|---|
| Generic case (Fig. 6) | < $1\times 10^{-4}$ from 15 keV; **0.5% at 1 keV** | — |
| $^{241}$Pu (7s valence) | Max ~1% in first 0.5 keV | Drops below $10^{-4}$ after 3 keV |

Effect of wave function change in lattice vs gas studied by Kolos et al. (1988) for tritium — no significant change found. Due to conservative error bar, **completely absorbed and can be neglected**.

## Calculator Implementation

> [!info] Not planned
> Chemical effects are NOT part of the calculator scope. They represent corrections on top of atomic calculations that only matter in specific experimental contexts (molecular tritium sources). The baseline calculator assumes isolated atoms.

For molecular experiments requiring these corrections:
- Use Saenz & Froelich formalism (§F.3) for screening modification
- Apply conservative error bars for exchange effects (~$10^{-4}$ from 15 keV onward)
- Q value shifts should be handled externally (experiment-specific input)

## Related Notes

- [[12-atomic-overlap]] — Atomic mismatch / Bahcall correction
- [[04-screening-correction]] — Baseline atomic screening $S(Z,W)$
