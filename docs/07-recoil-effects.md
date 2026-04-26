---
title: "07 – Recoil Effects"
date: 2026-04-26
tags:
  - beta-spectrum/correction/recoil
status: todo
aliases: [Nuclear Recoil, R_N, Weak Magnetism]
cssclasses: []
related_components: [[02-phase-space]], [[10-nuclear-structure]]
---

# Nuclear Recoil Effects — $R_N$, $Q$, and Induced Currents

## Status in Calculator

> [!todo] Not yet implemented
> Module: `beta_spectrum/components/recoil.py` exists but is **empty** (stub only). This correction needs to be added.

## Why Recoil Matters

The finite nuclear mass turns beta decay from a two-body into a three-body process, introducing recoil corrections that become significant at the $10^{-4}$ precision level. Additionally:
- **Weak magnetism** $b_{wm}$ must be included as experimental precision increases — otherwise it limits BSM sensitivity
- The oft-neglected **induced pseudoscalar interaction** may be significant

## Kinematic Recoil Correction $R_N(W,W_0,M)$ (Eq. 41–43)

$$R_N = 1 + \frac{r_0}{W} + \frac{r_1}{W^2} + r_2 W + \frac{r_3}{W^2} \tag{41}$$

### For Fermi (vector) transitions:
$$r_0^V = W_0^2/(2M^2) - 11/(6M^2), \quad r_1^V = W_0/(3M)$$
$$r_2^V = 2/M - 4W_0^2/(3M), \quad r_3^V = 16/(3M^2)$$

### For Gamow-Teller (axial) transitions:
$$r_0^A = -2W_0/(3M) - W_0^2/(6M^2) - 77/(18M^2), \quad r_1^A = -2/(3M) + 7W_0/(9M^2)$$
$$r_2^A = 10/(3M) - 28W_0/(9M^2), \quad r_3^A = 88/(9M^2)$$

### Mixed transitions:
Terms are weighted by $1/(1+\rho^2)$ and $1/(1+\rho^{-2})$ where $\rho = \mathcal{M}_{GT}/\mathcal{M}_F$.

> [!tip] Magnitude
> Order $10^{-5}$ to $10^{-3}$. See Wilkinson (1990) Fig. 1 for magnitude plot.

## Electromagnetic Recoil Correction $Q(Z,W,M)$ (Eq. 45)

The recoiling nucleus's Coulomb field differs from a static one:

$$Q \simeq 1 \mp \frac{\pi\alpha Z}{M}\left(1 + a\right)\frac{W_0 - W}{3W} \tag{45}$$

where $a = (1-\rho^2/3)/(1+\rho^2)$ is the β-ν correlation coefficient ($a=1$ for Fermi, $a=-1/3$ for GT).

> [!note] Typically negligible
> The effect amounts to at most a few percent of the typical error in phase space factor $f$ due to Q-value uncertainty — usually not limiting.

## Weak Magnetism and Induced Currents (Section VI)

The weak vertex receives QCD-induced form factors beyond pure V-A:

### Behrens-Bühring Form (Eq. 82):
$$J_\mu^{BB} = i\langle \bar{u}_p | C_V\gamma_\mu - f_M\sigma_{\mu\nu}q^\nu + if_S q_\mu $$
$$- \frac{C_A}{C_V}\left(\gamma_\mu\gamma^5 - f_T\sigma_{\mu\nu}\gamma^5 q^\nu\right) + i\frac{f_P}{C_V}\gamma^5 q_\mu | u_n \rangle$$

### Holstein Form (Eq. 83):
$$J_\mu^{HS} = i\langle \bar{u}_p | g_V\gamma_\mu - \frac{g_M-g_V}{2M}\sigma_{\mu\nu}q^\nu + i\frac{g_S}{2M}q_\mu $$
$$+ g_A\gamma^5\gamma_\mu - \frac{g_{II}}{2M}\sigma_{\mu\nu}\gamma^5 q^\nu + i\frac{g_P}{2M}\gamma^5 q_\mu | u_n \rangle$$

| Form Factor | Name | CVC/PCAC Status |
|---|---|---|
| $C_V = g_V$ | Vector (Fermi) | Conserved — from EM data |
| $f_M/C_V$ | Weak magnetism | **Determined by CVC** from M1 transitions |
| $f_S/g_S$ | Induced scalar | Excluded by CVC for Fermi; Second-class current |
| $C_A = g_A$ | Axial-vector (GT) | Measured in neutron decay ($g_A \approx 1.276$) |
| $f_T$ | Tensor | First-class if exists |
| $g_{II}$ | Induced tensor | Second-class current → experimental null results |
| $f_P/g_P$ | **Induced pseudoscalar** | PCAC prediction: $g_P(0) \approx -229$, but quenched ~80% in nuclei |

### CVC Constraint on Weak Magnetism (Eq. 87–91)

For mirror transitions and isovector M1 decays, weak magnetism is determined from electromagnetic data — no free parameters:
$$\mu_p - \mu_n = g_M(0)/\mu_N \quad (\text{from CVC})$$

### PCAC / Goldberger-Treiman (Eq. 92)

$$g_P(0) \approx -229 \quad (\text{free nucleon})$$

With ~5% higher-order corrections. In nuclear medium, $g_P$ is **quenched by up to 80%** — great care needed when evaluating pseudoscalar currents for nuclear decays.

## Spectral Shape Modification (Eq. 3)

Weak magnetism and Fierz term modify the spectrum shape:

$$N(W)\,dW \propto p W (W_0 - W)^2 \left[1 + \frac{\gamma m_e}{W} b_{\text{Fierz}} \pm \frac{4W}{3M} b_{wm}\right] dW \tag{3}$$

The weak magnetism term $\propto W/M$ grows with energy — opposite to Fierz term which falls as $1/W$. This distinct energy dependence allows experimental separation.

## Implementation Plan for Calculator

To implement recoil in `recoil.py`:
1. Create `RecoilCorrection` class inheriting from `SpectrumComponent`
2. Implement $R_N(W,W_0,M)$ per Eqs. (41)–(44) with transition-type parameter (Fermi/GT/mixed)
3. Optionally implement $Q(Z,W,M)$ for electromagnetic recoil
4. Accept weak magnetism form factor from nuclear data or CVC calculation
