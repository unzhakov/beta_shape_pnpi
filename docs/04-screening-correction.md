---
title: "04 – Screening Correction"
date: 2026-04-26
tags:
  - beta-spectrum/correction/screening
status: active
aliases: [Atomic Screening, S(Z,W)]
cssclasses: []
related_components: [[12-atomic-overlap]], [[13-chemical-effects]]
---

# Atomic Electron Screening Correction $S(Z,W)$

## Physical Origin

The Fermi function ([[01-fermi-function]]) assumes the β particle moves in the Coulomb field of a bare nucleus. In reality, atomic electrons partially screen this charge at large distances from the nucleus, modifying the effective potential experienced by the emitted electron/positron. This is particularly important at low kinetic energies where the β particle spends more time far from the nucleus.

## The Screening Correction

Screening is implemented as a ratio method:

$$S(Z,W) = \frac{F_{\text{screened}}(Z, W)}{F_0(Z, W - V_0/m_ec^2)} \tag{144}$$

where $V_0$ is the screening potential (binding energy of the outermost atomic electrons). The denominator evaluates the Fermi function at a shifted energy representing the unscreened case.

### Screening Potential

The screening potential $V_0$ depends on:
- Atomic number $Z$ (parent atom)
- Chemical/physical state (neutral atom, ionized molecule)
- Electronic configuration

For neutral atoms, $V_0$ is typically a few eV to tens of eV — small compared to nuclear potentials but significant at β energies below ~10 keV.

### Molecular Screening Correction

In molecular environments, the screening potential changes:

$$\Delta S_{\text{Mol}} \quad (\text{Eq. 176})$$

Since screening is typically at sub-percent level maximum, molecular deviations are expected at the **few $10^{-4}$** level — relevant for precision work down to 1 keV.

## Energy Dependence

- **Dominant at low energies** (T < ~50 keV): where the β particle's wave function is most affected by atomic electron cloud
- **Negligible at high energies**: $S(Z,W) \to 1$ as $W \gg V_0$
- The effect reverses sign between β⁻ and β⁺ (electrons attracted to/repelled from electron cloud differently)

## Shake-Up Effects on Screening

Atomic excitations during β decay modify screening:

> [!warning] Tritium case
> In tritium β⁻ decay, the screening potential can change by as much as **20%** because only a single electron contributes to atomic Coulomb interaction. For heavier systems, relative changes decrease but scale with $Z$. Neglecting shake-up introduces errors of order $\sim 1 \times 10^{-4}$.

Shake-off (continuum ejection) has smaller effect:
- With $\Delta Z_{\text{eff}} \approx 0.3$–$0.4$, yields ~**0.1%** relative change in screening potential → few $10^{-5}$ level in spectral shape
- **Negligible at target precision** for most nuclei

## Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/screening.py` → `ScreeningCorrection(FermiFunction)` class
> 
> Uses ratio method with smooth logistic switching function. The screening potential $V_0$ is passed as a parameter and depends on the parent atom's atomic number.

### API Usage

```python
from beta_spectrum import ScreeningCorrection, FermiFunction

screen = ScreeningCorrection(
    FermiFunction(Z_parent=90),  # bare nucleus reference
    V0=0.015  # screening potential in m_e units (~7.6 keV)
)
values = screen(W_grid)
```

## Size Comparison with Other Corrections

From Hayen et al. Table (Section VIII overview):

| Correction | Magnitude | Energy Range Most Relevant |
|---|---|---|
| Screening $S$ | $10^{-3}$ – $10^{-4}$ | T < 50 keV |
| Exchange $X$ | up to ~20% at lowest E | T < 10 keV |
| Molecular screening $\Delta S_{\text{Mol}}$ | few $10^{-4}$ | T < 5 keV (molecular) |

## Related Notes

- [[12-atomic-overlap]] — Atomic mismatch / overlap correction $r(Z,W)$, Bahcall correction
- [[13-chemical-effects]] — Molecular environment influence on screening and exchange
- [[05-exchange-correction]] — Electron exchange with atomic orbitals (often combined with screening)
