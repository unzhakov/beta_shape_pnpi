# Hayen et al. (2017) — High Precision Analytical Description of the Allowed β Spectrum Shape: Summary of Chunks 1–6

## Title, Authors, Abstract

**Title:** *High precision analytical description of the allowed β spectrum shape*

**Authors:** Leendert Hayen and Nathal Severijns (KU Leuven, Belgium); Kazimierz Bodek and Dagmara Rozpedzik (Jagiellonian University, Poland); Xavier Mougeot (CEA LIST, France)

**arXiv:** 1709.07530v2 [nucl-th], dated September 26, 2017

### Abstract Summary

The paper presents a fully analytical description of the allowed beta decay spectrum shape for ongoing and planned measurements aimed at searching for physics beyond the Standard Electroweak Model (BSM) and investigating weak magnetism recoil terms. Contributions from finite size corrections, mass effects, and radiative corrections are reviewed. A particular focus is placed on **atomic and chemical effects**, extending the existing description analytically. The effects of QCD-induced recoil terms are discussed with cross-checks across different theoretical formalisms. Corrections are derived for both Fermi and Gamow-Teller transitions. Calculated $f$ values agree with the most precise numerical results within the aimed-for precision. The paper stresses the need for accurate evaluation of weak magnetism contributions and notes the possible significance of the oft-neglected induced pseudoscalar interaction. Together with improved atomic corrections, an analytical description is presented **accurate to a few parts in $10^{-4}$ down to 1 keV** for low-to-medium $Z$ nuclei — extending previous work by nearly an order of magnitude.

---

## I. INTRODUCTION

### Motivation: Beta Decay as a Probe for New Physics

Beta decay played a pivotal role in uncovering the left-handed V-A weak interaction and the electroweak sector of the Standard Model (SM). In the LHC era, competitive BSM results are extracted from low-energy beta decay measurements through Effective Field Theories (EFT), where new physics at scale $\Lambda_{\text{BSM}} \gg \Lambda_{\text{LHC}}$ is parameterized.

The **Lee-Yang Hamiltonian** is generalized to include higher-dimensional SM field operator combinations. Searches for non-V-A manifestations focus on:
- Exotic scalar/tensor couplings (right-handed currents involving new heavy particles)
- Lorentz invariance violation via sidereal variations of observables

### The Fierz Interference Term

Scalar and tensor BSM coupling constants appear in the beta spectrum shape through the **Fierz interference term** (Jackson et al., 1957):

$$b_{\text{Fierz}} \simeq \pm \frac{1}{1 + \rho^2} \left[ \text{Re}\left(\frac{C_S + C'_S}{C_V}\right) + \rho^2 \, \text{Re}\left(\frac{C_T + C'_T}{C_A}\right) \right], \tag{1}$$

where $\rho = \mathcal{M}_{GT} / \mathcal{M}_F$ is the ratio of Gamow-Teller to Fermi matrix elements, and $C_{V,A}$ are SM vector/axial-vector coupling constants while $C_S, C_T$ are possible scalar/tensor couplings. For pure transitions:

$$b_{\text{Fierz},F} \simeq \pm \text{Re}\left(\frac{C_S + C'_S}{C_V}\right), \tag{2a}$$
$$b_{\text{Fierz},GT} \simeq \pm \text{Re}\left(\frac{C_T + C'_T}{C_A}\right). \tag{2b}$$

A precision of typically **0.5%** is required when determining $b_{\text{Fierz}}$ to improve limits on scalar/tensor coupling constants.

### Weak Magnetism and Spectral Shape Modification

QCD-induced form factors — dominated by the **weak magnetism term** $b_{wm}$ — must be included as experimental precision increases, since they otherwise limit BSM sensitivity:

$$N(W)\,dW \propto pW(W_0 - W)^2 \left[ 1 + \frac{\gamma m_e}{W} b_{\text{Fierz}} \pm \frac{4W}{3M} b_{wm} \right] dW, \tag{3}$$

where:
- $p$, $W$, $W_0$ are the beta particle momentum, total energy, and endpoint energy respectively
- $\gamma = \sqrt{1 - (\alpha Z)^2}$ with $\alpha$ the fine structure constant and $Z$ the daughter atomic number
- $m_e$ is the electron rest mass and $M$ is the average mother/daughter nucleus mass

To search for new physics or establish weak magnetism effects, a beta spectrum description reliable at the **$10^{-4}$ level** is required.

### Previous Theoretical Approaches

Three major formalisms are reviewed:
1. **Behrens-Bühring** (and collaborators) — rigorous but relies on numerical solutions of the Dirac equation
2. **Holstein** — analytical approach with clean notation
3. **Wilkinson** — analytical parametrizations of dominant effects

The paper's goal is to combine the computational accuracy of Behrens-Bühring with the transparency of Holstein, achieving $10^{-4}$ precision analytically.

### Atomic and Chemical Effects (New Contributions)

Special attention improves the analytical description of:
- **Screening effects** — combining Bühring's theoretical work with high-precision atomic potentials for all $Z$ values; an analytical fit is proposed
- **Atomic exchange effect** — based on state-of-the-art numerical calculations, its contribution can exceed 20% in low-energy regions for higher $Z$
- **Shake-up and shake-off processes** — reviewed for their influence on atomic effects
- **Molecular effects** — explored for the first time at this precision level

---

## II. COMPLETE EXPRESSION

The fully differential allowed beta spectrum shape, including all corrections needed for $10^{-4}$ precision, is:

$$N(W)\,dW = \frac{G_V^2 V_{ud}^2}{2\pi^3} F_0(Z,W) L_0(Z,W) U(Z,W) DFS(Z,W,\beta^2) R(W,W_0) R_N(W,W_0,M)$$
$$\times Q(Z,W) S(Z,W) X(Z,W) r(Z,W) C(Z,W) D_C(Z,W,\beta^2) \, pW(W_0 - W)^2 dW \tag{4}$$

$$\equiv \frac{G_V^2 V_{ud}^2}{2\pi^3} K(Z,W,W_0,M) \, A(Z,W) \, C'(Z,W) \, pW(W_0-W)^2 dW.$$

### Definitions of All Factors

| Symbol | Name | Description |
|---|---|---|
| $F_0(Z,W)$ | Point charge Fermi function | Coulomb interaction between beta particle and daughter nucleus (Eq. 6) |
| $L_0(Z,W)$ | Electrostatic finite size correction | Corrections to Fermi function when evaluating at origin for uniform sphere |
| $U(Z,W)$ | More realistic charge distribution correction | Difference between actual nuclear shape and uniform sphere |
| $DFS(Z,W,\beta^2)$ | Deformation correction | Effect of axially symmetric nuclear deformation on the Fermi function |
| $R(W,W_0)$ | Outer radiative corrections | QED corrections (Sec. V) |
| $R_N(W,W_0,M)$ | Finite nuclear mass recoil correction | Kinematic effects from finite nucleus mass (Sec. IV.B.1) |
| $Q(Z,W,M)$ | Recoil electromagnetic correction | Coulomb field modification due to nuclear recoil (Sec. IV.B.2) |
| $S(Z,W)$ | Screening correction | Atomic electron screening effects (Sec. VII.A) |
| $X(Z,W)$ | Atomic exchange effect | Electron exchange with atomic orbitals (Sec. VII.B) |
| $r(Z,W)$ | Atomic mismatch / overlap | Alternative atomic excitation correction (Sec. VII.D) |
| $A(Z,W)$ | Combined atomic factors | $S \cdot X \cdot r$ |
| $C(Z,W)$ | Nuclear structure sensitive correction | Contains all nuclear form factors and matrix elements (Sec. VI) |
| $D_C(Z,W,\beta^2)$ | Nuclear deformation correction | Nuclear shape effects in the convolution terms (Sec. VI.F) |

### Key Conventions

- **Natural units**: $c = \hbar = m_e = 1$ throughout all formulae
- **Z convention**: $Z$ is always positive; $\beta^-$ and $\beta^+$ are distinguished explicitly via upper/lower signs ($\mp$)
- $W = E/m_ec^2 + 1$ (total beta energy in units of electron rest mass)
- $p = \sqrt{W^2 - 1}$ (momentum in units of $m_ec$)
- $G_V$: vector coupling strength in nuclei; $V_{ud} = \cos^2\theta_C$

### Terminology Clarification: "Finite Size Effects"

The authors distinguish two types that different authors conflate under "finite size":

1. **Electrostatic finite size** ($L_0, U, DFS$): Stems from the electric potential difference when moving from point charge to realistic nuclear shape. These are corrections applied *after* evaluating the electron wave function at the origin.

2. **Convolution finite size**: Results from integrating leptonic wave functions over the nuclear volume. This is sensitive to nuclear structure and appears in $C(Z,W)$. The authors note that different communities use "finite size" to refer to entirely different corrections, and aim to clarify this distinction.

---

## III. FERMI FUNCTION

### The Point Charge Fermi Function

The transition matrix element (to first order in the weak Hamiltonian with all-order Coulomb interactions resummed) is:

$$\mathcal{M}_{fi} = -2\pi i \delta(E_f - E_i) \langle f | T \exp\left(-i\int_0^\infty dt H_Z(t)\right) $$
$$\times H_\beta(0) T \exp\left(-i\int_{-\infty}^0 dt H_{Z'}(t)\right) | i \rangle, \tag{5}$$

The typical approximation is to keep only the final-state Coulomb interaction (daughter nucleus charge $Z$), with corrections for radiative effects and atomic influences.

**Fermi's point charge Fermi function** (Eq. 6):

$$F_0(Z,W) = 4(2pR)^{2(\gamma-1)} e^{\pi y} \frac{|\Gamma(\gamma + iy)|^2}{(\Gamma(1+2\gamma))^2}, \tag{6}$$

with:
$$\gamma = \sqrt{1 - (\alpha Z)^2}, \quad y = \pm \alpha Z W / p. \tag{7}$$

Here $R$ is the cutoff radius (representing daughter nucleus radius) necessitated by divergence at origin for point charge, though it has no real physical significance — most of its dependence is absorbed into corrections.

### Alternative Fermi Function Definitions

The paper reviews multiple conventions:

- **Wilkinson / Blatt & Weisskopf / Fano / Konopinski style**: Uses prefactor $(1+\gamma)/2$ instead of 4, absorbing it into $L_0$.
- **Behrens-Bühring / Schopper style** (followed by this paper): Keeps Eq. (6) with factor 4, and includes all additional corrections in $L_0$ and subsequent factors.

### Finite Size Correction from Power Expansion

Expanding the radial Dirac wave functions near the origin:

$$\left\{\frac{f_\kappa(r)}{g_\kappa(r)}\right\} = \alpha_\kappa \{(2|\kappa|-1)!!\}^{-1}(pr)^{|\kappa|-1} \sum_{n=0}^\infty \frac{a_{\kappa n}}{b_{\kappa n}} r^n, \tag{12}$$

the Coulomb information is encoded in **Coulomb amplitudes** $\alpha_\kappa$. The dominant electron component at the origin becomes:

$$F_0 L_0 = \frac{\alpha_{-1} + \alpha_1^2}{2p^2}. \tag{13}$$

This is valid regardless of nuclear charge distribution shape.

### Behrens-Jänecke Finite Size Correction

An alternative Fermi function that includes finite size:

$$F_{BJ}(Z,W) = F_0(Z,W)\left[1 \mp \frac{\alpha Z W R}{15} + ...\right], \tag{10}$$

The additional term traces to the expansion of $f_1$ and $g_{-1}$ inside the nucleus (Huffaker and Laird, 1967):

$$\left\{\frac{g_{-1}(r)}{f_1(r)}\right\}^2 \propto [p F_0(Z,W)]^{1/2}\left[1 - \frac{13}{30}\left(\frac{1}{2R}\right)^2 + \alpha Z W R - (pr)^2\right], \tag{11}$$

interpreted as a rough electrostatic finite size correction. The error from neglecting the $(1+\gamma)/2$ prefactor is $\sim 0.5\%$ for $Z=20$.

---

## IV. FINITE MASS AND ELECTROSTATIC FINITE SIZE EFFECTS

### A. Electrostatic Finite Nuclear Size Corrections

#### 1. L₀(Z, W) — Uniformly Charged Sphere

When the nucleus has finite size rather than being a point charge, electron/positron wave functions become finite at the center. The correction $L_0$ is introduced for use with the analytical Fermi function in Eq. (6).

The nuclear radius is adjusted to match experimental $\langle r^2 \rangle^{1/2}$:

$$R = \sqrt{5/3} \, \langle r^2 \rangle^{1/2}. \tag{14}$$

**Low-Z expansion** (Behrens and Bühring):

$$L_0 \simeq \frac{1+\gamma}{2}\left[1 \mp \alpha Z W R + (\alpha Z)^2 - \gamma \frac{\alpha Z R}{2W}\right] $$
$$\simeq 1 \mp \alpha Z W R + (\alpha Z)^2 - \frac{13}{60} - \frac{\alpha Z R}{2W}. \tag{15}$$

**Wilkinson's analytical parametrization** (accurate to $10^{-4}$ for $p \leq 45$, $|Z| \leq 60$):

$$L_0(Z,W) = 1 + \frac{13}{60}(\alpha Z)^2 \mp \frac{\alpha Z W R(41 - 26\gamma)}{[15(2\gamma-1)]} $$
$$\mp \frac{\alpha Z R \gamma(17 - 2\gamma)}{[30W(2\gamma-1)]} + \frac{R}{W}\sum_{n=0}^5 a_n (WR)^n + 0.41(R-0.0164)(\alpha Z)^{4.5}, \tag{16}$$

with coefficients expanded as:
$$a_n = \sum_{x=1}^6 b_{x,n} (\alpha Z)^x. \tag{17}$$

Coefficients $b_{x,n}$ are tabulated in **Tables I (electrons)** and **II (positrons)**. For positrons, signs of odd powers of $Z$ are flipped. The effect ranges from a few 0.1% up to several percent — everywhere highly significant.

#### 2. U(Z, W) — More Realistic Charge Distributions

The uniformly charged sphere is an approximation; actual nuclear charge distributions include:
- **Modified Gaussian**: $\rho_{MG}(r) = N_0\left[1 + A\left(\frac{r}{a}\right)^2\right]e^{-(r/a)^2}$ (Eq. 18), with normalization $N_0 = \frac{8 a^{-3} \pi^{-1/2}}{2+3A}$ and $a = R\sqrt{\frac{5}{2}(2+5A)/(2+3A)}$ (Eqs. 19–20)
- **Fermi distribution**: $\rho_F(r) \propto \left[1 + \exp((r-c)/a)\right]^{-1}$ with $a \simeq 0.55$ fm

The correction from replacing uniform sphere with a more realistic distribution is:

$$\frac{U(Z,W)}{L_0} = \frac{\alpha_{02}' + \alpha_{-1}'^2}{\alpha_{1}^2 + \alpha_{-1}^2} = \frac{B_{10^{-2}}' + B_{-1}'^2}{B_{1^{-2}} + B_{-1}^2}, \tag{23}$$

neglecting $\mathcal{O}(WR^2)$ terms, where $B_{\pm 1}$ involves wave functions inside the nucleus (Eq. 24).

**Analytical approximation for low Z** (Eq. 25):
$$U(Z,W) \approx 1 + \alpha Z W R \Delta_1 + \frac{\gamma}{W} \alpha Z R \Delta_2 + (\alpha Z)^2 \Delta_3 - (WR)^2 \Delta_4, \tag{25}$$

where $\Delta_i$ are linear combinations of potential coefficient differences $\Delta v_n = v_n^0 - v_n'$ (Eqs. 26a–d).

**Wilkinson's Fermi distribution fit**:
$$U(Z,W) = 1 + \sum_{n=0}^\infty a_n p^n, \tag{29}$$
with coefficients:
$$a_0 = -5.6\times 10^{-5} \mp 4.94\times 10^{-5} Z + 6.23\times 10^{-8} Z^2,$$
$$a_1 = 5.17\times 10^{-6} \pm 2.517\times 10^{-6} Z + 2.00\times 10^{-8} Z^2,$$
$$a_2 = -9.17\times 10^{-8} \pm 5.53\times 10^{-9} Z + 1.25\times 10^{-10} Z^2, \tag{30}$$

The effect of $U(Z,W)$ is $\sim 0.1\%$ for medium-high masses — cannot be neglected. The three-parameter Fermi ("wine bottle") distribution can use Eq. (29) as first approximation with Eq. (25) describing differences between two Fermi distributions.

#### 3. DFS(Z, W, β²) — Nuclear Deformation

For axially symmetric deformations:
$$R(\theta,\phi) = R_0\left[1 + \beta_2 Y_{20}(\theta,\phi)\right], \tag{32}$$

with $\beta_2 > 0$ prolate, $\beta_2 < 0$ oblate. The intrinsic electric quadrupole moment:
$$Q_0 = \frac{\alpha}{\pi}\sqrt{\frac{5}{3}} R Z \beta_2(1 + 0.16\beta_2). \tag{33}$$

The deformation correction to the Fermi function is defined as a ratio of angle-averaged corrections:
$$DFS(Z,W,\beta) = \frac{L_0^*(Z,W)}{L_0(Z,W)}, \tag{40}$$

where $L_0^*$ involves integrating over a continuous superposition of uniformly charged spheres (Eq. 39). The effect can reach several parts in $10^3$ but is partially canceled by the nuclear structure deformation correction $D_C$. **Figure 1** shows DFS for $\beta_2 = 0.2$, revealing an energy dependence reversal between $\beta^-$ and $\beta^+$ decay.

### B. Finite Nuclear Mass Corrections

#### 1. Rₙ(W, W₀, M) — Recoil Kinematic Correction

The finite nuclear mass turns beta decay from two-body to three-body process:
$$R_N(W,W_0,M) = 1 + r_0 + r_1/W + r_2 W + r_3/W^2. \tag{41}$$

**For Fermi (vector) transitions**:
$$r_0^V = W_0^2/(2M^2) - 11/(6M^2), \tag{42a}$$
$$r_1^V = W_0/(3M), \tag{42b}$$
$$r_2^V = 2/M - 4W_0^2/(3M), \tag{42c}$$
$$r_3^V = 16/(3M^2). \tag{42d}$$

**For Gamow-Teller (axial) transitions**:
$$r_0^A = -2W_0/(3M) - W_0^2/(6M^2) - 77/(18M^2), \tag{43a}$$
$$r_1^A = -2/(3M) + 7W_0/(9M^2), \tag{43b}$$
$$r_2^A = 10/(3M) - 28W_0/(9M^2), \tag{43c}$$
$$r_3^A = 88/(9M^2). \tag{43d}$$

For mixed transitions, terms are weighted by $1/(1+\rho^2)$ and $1/(1+\rho^{-2})$ (Eq. 44). This is a small effect of order $10^{-5}$ to $10^{-3}$. See Wilkinson (1990) Fig. 1 for magnitude plot.

#### 2. Q(Z, W, M) — Recoil Electromagnetic Correction

The recoiling nucleus's Coulomb field differs from a static one:
$$Q(Z,W,M) \simeq 1 \mp \frac{\pi \alpha Z}{M}\left(1 + a\right)\frac{W_0 - W}{3W}, \tag{45}$$

where $a = (1-\rho^2/3)/(1+\rho^2)$ is the $\beta$-$\nu$ correlation coefficient ($a=1$ for Fermi, $a=-1/3$ for Gamow-Teller). The effect amounts to at most a few percent of the typical error in phase space factor $f$ due to Q-value uncertainty — typically negligible.

---

## V. RADIATIVE CORRECTIONS

### Overview: Inner vs. Outer Corrections

The radiatively corrected spectrum is:
$$\frac{d\Gamma}{dW} = \frac{d\Gamma_0}{dW}(1 + \Delta_R^{V/A})(1 + \delta_R(W,W_0)) \tag{46}$$

- **Inner corrections** $\Delta_R$: Nuclear-independent, absorbed into effective coupling constants. For vector case: $\Delta_V^R = (2.361 \pm 0.038)\%$ (Marciano and Sirlin, 2006).
- **Outer corrections** $\delta_R(W,W_0)$: Energy-dependent, nucleus-dependent — the focus of this section, contributing a few percent.

### A. Total Spectral Influence (Outer Radiative Corrections)

$$R(W,W_0) = 1 + \delta_R(W,W_0). \tag{47}$$

The correction is expressed as:
$$\delta_R = \frac{\alpha}{2\pi}\left[g(W_0,W) - 3\ln\left(\frac{m_p}{2W_0}\right)\right] + Z\alpha^2\sum_i \Delta_i(W) + Z^2\alpha^3 \delta_{3h} + ..., \tag{48}$$

where:
- $L(2W_0, m_p) = 1.026725\left[1 - \frac{2\alpha}{3\pi}\ln(2W_0)\right]$ (Eq. 49) — included in the constant term normalization

#### Order α Correction — Sirlin's g(W₀, W) Function:

$$g(W_0,W) = 3\ln(m_p) - \frac{3}{4} + \frac{4\beta}{1+\beta}\left[4\tanh^{-1}\beta - 1\right] $$
$$+ \frac{W_0-W}{\beta}\left(-\frac{3}{3W} + \ln[2(W_0-W)]\right) $$
$$- \frac{(W_0-W)^2}{6W^2} + 4\frac{\tanh^{-1}\beta}{\beta}\cdot 2(1+\beta^2) - \frac{(W_0-W)^2}{6W^2} - 4\tanh^{-1}\beta, \tag{50}$$

where $\beta = p/W$ and:
$$L_s(x) = \int_0^x \frac{\ln(1-t)}{t} dt = -\sum_{k=1}^\infty \frac{x^k}{k^2} \equiv -Li_2(x), \tag{51}$$

the **Spence function (dilogarithm)**. Its large $W_0$ limit:
$$g(W_0 \to \infty, W) = 3\ln\left(\frac{M_N}{2W_0}\right) + \frac{81}{10} - \frac{4\pi^2}{3}. \tag{52}$$

Equation (50) is universal — same for electrons and positrons, Fermi and Gamow-Teller, independent of nucleus. It has a **logarithmic divergence at $W = W_0$** related to soft real photon emission.

#### Soft Photon Resummation:
To handle the endpoint divergence, Eq. (53) is used:
$$t(\beta)\ln(W_0 - W) \to (W_0-W)^{t(\beta)} - 1, \tag{53}$$
where:
$$t(\beta) = \frac{2\alpha}{\pi}\left[\frac{\tanh^{-1}\beta}{\beta} - 1\right]. \tag{54}$$

The change in tritium spectrum is negligible due to low endpoint, but for higher energies the $f t$ value shift can reach several parts in $10^4$.

#### Order Zα² Correction:
Three dominant Feynman diagrams (Fig. 3): vacuum polarization of bremsstrahlung photon and two electron propagator renormalizations. The axial vector component contributes specifically to $\Delta_3$.

$$\delta_2(Z,W) = Z\alpha^2 \sum_{i=1}^4 \Delta_i(W). \tag{55}$$

The leading terms $\Delta_1$ depend on nuclear form factor $F(q^2)$ and are split into nucleus-independent ($\Delta_1^0$) and nucleus-dependent ($\Delta_1^F$) parts:
- Non-relativistic approximation (Eq. 56): energy-dependent but nucleus-independent terms combined with $\Delta_4$
- For uniformly charged sphere, Eqs. (57)–(62) give explicit integral forms for $\Delta_1^F$, $\Delta_2$, and $\Delta_3$

For modified Gaussian model results are similar; differences between models for superallowed Fermi transitions ($^{10}$O to $^{54}$Co) are $< 10^{-4}$.

#### Order Z²α³ Correction — Sirlin's Heuristic Formula:
$$\delta_{3h} = Z^2\alpha^3 \left[a\ln\left(\frac{\Lambda}{W}\right) + bf(W) + \frac{4\pi}{3}g(W) - 0.649\ln(2W_0)\right], \tag{63}$$

with $a = \frac{\pi}{3} - \frac{\sqrt{3}}{2}$ (Eq. 64), $b = \frac{4}{3\pi^4}(\pi^2 - \gamma_E) - \frac{\pi^2}{18}$ (Eq. 65), $f(W) = \ln(2W)-5/6$ (Eq. 66), $g(W)$ in Eq. (67).

#### Higher Order Corrections — Wilkinson Estimate:
$$\delta_{Z^n\alpha^{n+1}} \approx Z^n\alpha^{n+1} K_{nm} \ln^{m-n}\left(\frac{\Lambda}{W}\right), \tag{68}$$

with $K_{nm} \approx 0.5$ on average, giving summed contribution:
$$\delta^{\text{higher}} = \sum_{n=3}^\infty Z^n\alpha^{n+1} = \frac{Z^3\alpha^4}{1-Z\alpha}. \tag{69}$$

These become relevant for higher $Z$.

---

## B. Neutrino Radiative Corrections

Even though the neutrino has no direct electromagnetic interaction, it is indirectly influenced via virtual photon exchange and energy conservation from inner bremsstrahlung (relevant for reactor neutrino oscillation studies). The total $\mathcal{O}(\alpha)$ radiative correction to the antineutrino spectrum:

$$R_\nu(Z,W,W_0) = 1 + \frac{\alpha}{2\pi}\left[3\ln\left(\frac{m_p}{m_e}\right) + \frac{23}{4} - Li_2\left(-\frac{8}{\hat{\beta}}\right)\right.$$
$$+ 8\left(\frac{2\hat{W}\hat{\beta}}{\hat{\beta}-1}\right)^{-1} + 4\frac{\tanh^{-1}\hat{\beta}}{\hat{\beta}}\frac{7+3\hat{\beta}^2}{8} - 2\tanh^{-1}\hat{\beta}) \tag{70}$$

where $\hat{W} = W_0 - W_\nu$, $\hat{p} = \sqrt{\hat{W}^2 - m_e^2}$, $\hat{\beta} = \hat{p}/\hat{W}$ (Eq. 71). The correction is much smaller than for the electron spectrum.

---

## C. Radiative Beta Decay — Internal Bremsstrahlung

An additional photon in the final state:
$$n \to p + e^- + \bar{\nu}_e + \gamma, \tag{72}$$

Understood classically as radiation from sudden beta particle acceleration (Bloch 1936; Knipp and Uhlenbeck). The photon ejection probability per unit energy:

$$\Phi(W,\omega) = \frac{\alpha p}{\pi\omega p'}\left[\frac{W+W'^2}{p^2}\ln\left(\frac{W+p}{W'-p}\right) - 2\right], \tag{73}$$

with $W' = W_0-\omega$ and $p'=\sqrt{W'^2-1}$ the redefined energy/momentum. The photon spectrum affecting the beta spectrum is found by averaging over intermediate states (Eq. 74):
$$S(\omega) = \int_{1+\omega}^{W_0} N(W')\Phi(W',\omega)dW'.$$

Detection of internal bremsstrahlung reveals underlying weak interaction physics through correlations with lepton momenta and circular polarization. Experimental branching ratios: $BR_{\beta\gamma} \sim 3\times 10^{-3}$ for neutron, $\sim 2\times 10^{-3}$ for $^{32}$P.

---

## VI. NUCLEAR STRUCTURE EFFECTS — THE SHAPE FACTOR (Partial)

### Form Factors and Multipole Expansion

Nuclear structure effects are described via the elementary particle approach with form factors $F_K$ absorbing all nuclear structure information as functions of $q^2 = (p_f - p_i)^2$:

$$\langle f | V_\mu - A_\mu | i \rangle \propto \sum_{KMs} \sum_{L=K-1} (-1)^{J_f-M_f+M} (-i)^L F_{KLs}(q^2) $$
$$\times \sqrt{4\pi(2J_i+1)} \begin{pmatrix} J_f & K & J_i \\ -M_f & M & M_i \end{pmatrix} \frac{(qR)^L}{(2L+1)!!}. \tag{78}$$

Since $(qR)^2$ is very small for beta decay, the expansion typically stops after the first order. The formalism connects to both Behrens-Bühring and Holstein approaches, utilizing Behrens-Bühring's computational machinery with Holstein's clean notation.

---

# Hayen et al. (2017) — High Precision Analytical Description of the Allowed β Spectrum Shape: Summary of Chunks 7 (Section VI.A–C)

## A. Introduction to Nuclear Structure Effects

### The Beta Decay Hamiltonian and Current Decomposition

The beta decay Hamiltonian is constructed as a current-current interaction:

$$H_\beta(x) = \frac{G_V \cos\theta_C}{\sqrt{2}} [J^\mu(x)L_\mu(x) + \text{h.c.}], \tag{75}$$

where $\theta_C$ is the Cabibbo angle, $J^\mu(x)$ is the nuclear current, and $L^\mu(x)$ is the lepton current:

$$L^\mu(x) = i\bar{\varphi}_e(x)\gamma^\mu(1+\gamma^5)v_\nu(x). \tag{77}$$

Here $\varphi_e$ is the solution of the beta particle wave function in a Coulomb potential. The strong interaction requires generalizing the nuclear current from pure V-A:

$$J^{\mu\dagger}(x) = \langle f | V^\mu(x) + A^\mu(x) | i \rangle, \tag{76}$$

The transition matrix element (TME) is constrained by angular momentum coupling rules for $J_i^{\pi_i} \to J_f^{\pi_f}$ transitions, allowing decays of different multipole orders $K$ within the vector triangle $(J_i, J_f, K)$. The general shape factor formulation groups all non-phase-space/non-Fermi-function terms into:

$$N(W)\,dW \propto pW(W_0-W)^2 F(Z,W) C(Z,W). \tag{79}$$

**Key distinction:** $C(Z,W)$ encompasses two distinct effects:
1. **Spatial variation** of leptonic and nuclear wave functions inside the finite nuclear volume
2. **Raw nuclear structure** (form factors, matrix elements)

All corrections from previous sections (electrostatic finite size, radiative corrections, etc.) are *not* included in $C(Z,W)$ — they are treated as separate higher-order corrections.

### Two Formalisms Compared: Behrens-Bühring vs. Holstein

The paper discusses two major approaches to nuclear structure effects:
- **Behrens-Bühring (BB):** Couples leptonic and nuclear parts together; provides the most precise results through rigorous expansion in $r^2$, $(m_eR)^a$, $(WR)^b$, $(\alpha Z)^c$ introducing tabulated functions $I(k_e, m, n, \rho)$
- **Holstein (HS):** Decouples leptonic and nuclear parts; uses manifestly covariant form factors with cleaner symmetry properties

The authors adopt the BB computational approach but present final results in Holstein notation for clarity. See Appendix E for detailed comparison.

---

## B. Form Factors

### 1. Impulse Approximation

Form factors $F_{KLs}(q^2)$ describe nuclear structure model-independently as functions of the only Lorentz invariant $q^2 = (p_f - p_i)^2$. Using the **impulse approximation**, nucleons are treated as non-interacting particles coupling to the weak vertex as free particles:

$$O_{KLs} = \sum_{\alpha,\beta} \langle \alpha | O_{KLs} | \beta \rangle a^\dagger_\alpha a_\beta, \tag{80}$$

with nuclear matrix elements decomposing into:

$$\langle f | O_{KLs} | i \rangle = \sum_{\alpha,\beta} \langle \alpha | O_{KLs} | \beta \rangle \langle f | a^\dagger_\alpha a_\beta | i \rangle. \tag{81}$$

The one-body density matrix elements $\langle f | a^\dagger_\alpha a_\beta | i \rangle$ are calculable via shell model methods.

### 2. Induced Currents

QCD influences from within the nuclear potential add Lorentz-invariant terms to the weak vertex. In both BB and HS formalisms (completely equivalent):

**Behrens-Bühring form:**
$$J_\mu^{BB} = i\langle \bar{u}_p | C_V\gamma_\mu - f_M\sigma_{\mu\nu}q^\nu + if_S q_\mu $$
$$- \frac{C_A}{C_V}\left(\gamma_\mu\gamma^5 - f_T\sigma_{\mu\nu}\gamma^5 q^\nu\right) + i\frac{f_P}{C_V}\gamma^5 q_\mu | u_n \rangle, \tag{82}$$

**Holstein form:**
$$J_\mu^{HS} = i\langle \bar{u}_p | g_V\gamma_\mu - \frac{g_M-g_V}{2M}\sigma_{\mu\nu}q^\nu + i\frac{g_S}{2M}q_\mu $$
$$+ g_A\gamma^5\gamma_\mu - \frac{g_{II}}{2M}\sigma_{\mu\nu}\gamma^5 q^\nu + i\frac{g_P}{2M}\gamma^5 q_\mu | u_n \rangle, \tag{83}$$

where $\sigma_{\mu\nu} = -\frac{i}{2}[\gamma_\mu,\gamma_\nu]$ and all $C,g$ functions depend on $q^2$. The two approaches are equivalent through redefinition of form factors.

**Example — induced scalar modification (Eq. 84):**
$$^{(0)}_{F000} \to {}^{VM}_{000} \pm \frac{f_S}{R}\frac{6}{5}(W_0 R \pm \alpha Z)\,{}^{VM}_{000}, \tag{84}$$

### Table III: Nuclear Matrix Elements for Allowed Transitions

| Form Factor (BB) | Cartesian Form | $\Delta J$ | Forbidden | Type |
|---|---|---|---|---|
| $V_{F000}^{(0)}$ | $+g_V \mathbf{1}$ | 0 | — | **Allowed** |
| $A_{F101}^{(0)}$ | $\mp g_A \boldsymbol{\sigma}$ | 0,1 | $0^- \to 0^-$ | **Allowed** |
| $V_{F000}^{(1)}$ | $+g_V (\mathbf{r}/R)^2$ | 0 | — | Main correction |
| $A_{F101}^{(1)}$ | $\mp g_A \frac{\mathbf{r}}{R}\left(\frac{R}{\sqrt{2(M_N R)^2}} + \boldsymbol{\sigma}\right)$ | 0,1 | $0^- \to 0^-$ | Main correction |
| $A_{F121}^{(0)}$ | $\mp g_A\sqrt{\frac{3}{2}}\frac{\sqrt{5}}{R^2}[\sqrt{4\pi}\mathbf{Y}_2(\hat{\mathbf{r}})\otimes\boldsymbol{\sigma}]_1 \pm (2M_N R)^{-2}$ | 0,1 | $0^- \to 0^-$ | — |
| $V_{F011}^{(0)}$ | $+g_V i\frac{\mathbf{r}}{R}$ | 0 | — | **Relativistic** |
| $V_{F111}^{(0)}$ | $-g_V\sqrt{\frac{3}{2}}\boldsymbol{\alpha}\times\frac{\mathbf{r}}{R} - \sqrt{\frac{M_N}{2}}\frac{g_2^M}{NR}$ | 0,1 | $0^- \to 0^-$ | **Relativistic** |
| $A_{F110}^{(0)}$ | $\pm g_A\sqrt{3}\gamma^5\frac{\mathbf{r}}{R} \pm (2M_N R)^{-2}W_0 R [\cdots] \pm \frac{g_P}{\sqrt{56}}\alpha Z$ | 0,1 | $0^- \to 0^-$ | **Relativistic** |

### Table IV: Holstein Form Factors in Impulse Approximation

| Form Factor | Formula (Impulse Approx.) | Remark | Type |
|---|---|---|---|
| $a$ | $\sim g_V M_F$ | $g_V = 1$ (CVC, Ademollo-Gatto) | Vector |
| $e$ | $\approx g_V(M_F \pm A g_S)$ | $e=0$ (CVC, SCC) | Vector |
| $b$ | $\approx b = A(g_M M_{GT} + g_V M_L)$ | $g_M \approx 4.706$ | Vector |
| $c_1$ | $\approx c_1 = g_A M_h^{GT}$ | $g_A \to g_{A,\text{eff}} = 1$ (Towner, 1987) | Axial vector |
| $c_2$ | $\approx \frac{i}{6}g_A M_{\sigma r^2} + \frac{\sqrt{1}}{10}M_{ky}$ | $c_2 \sim R^2$ | Axial vector |
| $d$ | $\approx d = A(g_A M_{\sigma L} \pm g_{II} M_{GT})$ | $g_{II} \sim g_T \approx 0$ (SCC); $d_I=0$ for analog states | Axial vector |
| $h$ | $\approx -\sqrt{\frac{2}{10}}M_2^2 g_A M_{1y} - A^2 g_P M_{GT}$ | $g_{P,\text{free}} \approx -229 \to g_{P,\text{eff}} = ?$ | Axial vector |

### Table V: Nuclear Matrix Elements (Calaprice et al., 1977b)

| Matrix Element | Operator Form |
|---|---|
| $M_F$ | $\langle\beta\|\sum\tau_i^\pm\| \alpha_i\rangle$ |
| $M_{GT}$ | $\langle\beta\|\sum\tau_i^\pm \vec{\sigma}_i\| \alpha_i\rangle$ |
| $M_L$ | $\langle\beta\|\sum\tau_i^\pm \vec{l}_i\times\vec{\sigma}_i\| \alpha_i\rangle$ |
| $M_{\sigma r^2}$ | $\langle\beta\|\sum\tau_i^\pm \vec{\sigma}_i r_i^2\| \alpha_i\rangle$ |
| $M_{ky}^{(0)}$ | $\frac{1}{6}\sqrt{\frac{4\pi}{5}}\langle\beta\|\sum\tau_i^\pm r_i^2 C_1^{(0)}Y_{2n}(\hat{\mathbf{r}}_i)\| \alpha_i\rangle$ |

### 3. Validity of the Impulse Approximation

The impulse approximation neglects many-body effects (meson exchange, core polarization). Key findings:
- **Vector matrix element $V_{F000}$:** The impulse approximation is valid because the vector current is conserved (CVC)
- **Axial current:** Breakdown occurs due to PCAC and pion field coupling — corrections are significant

**Two categories of nuclear structure corrections:**
1. **Meson exchange effects** — long history (Blin-Stoyle et al., Chemtob & Rho, Delorme & Rho)
2. **"Nuclear" effects** — core polarization, relativistic effects, configuration mixing from insufficient wave function knowledge

These distinctions are artificial since the nuclear potential incorporates mesonic degrees of freedom (Wilkinson, 1974). For **allowed decays**, impulse approximation results largely hold due to cancellations between core-polarization and meson exchange (Morita, 1985), though individual effects can reach ~40%.

**Quenching of coupling constants:**
- $g_A$ is renormalized ad hoc: $g_{A,\text{eff}} = 1.1$ for sd-shell nuclei, unity for fp-shell (Martínez-Pinedo et al., 1996), and even $g_A \approx 0.9$ from highly forbidden decays
- The **induced pseudoscalar coupling** $g_P$ can be quenched by up to **80%** in nuclear matter — significantly more than the connected quenching of $g_A$
- Ab initio methods (Green's function Monte Carlo) circumvent quenching but fall outside scope

---

## C. Symmetries

### 1. Conserved Vector Current (CVC)

The CVC hypothesis (Feynman and Gell-Mann, 1958) provides powerful constraints:

**Temporal component:**
$$\langle J_f M_f | V^0(0) | J_i M_i \rangle = {}^{V}F_{000}(q^2)\delta_{J_i J_f}\delta_{M_i M_f} $$
$$= \sqrt{(T_i \pm T_{3i})(T_i \mp T_{3i} + 1)}\,\delta_{J_i J_f}\delta_{M_i M_f}, \tag{87-88}$$

This implies CVC **excludes induced scalar currents** ($e=0$ in Table IV).

**Relation between vector form factors (Eq. 89):**
$$-2N\,{}^{V}F_{011}^{N-1} = \left(W_0 \mp (m_n - m_p)\right)R\,{}^{V}F_{000}^N $$
$$\pm \alpha Z \int_R \!\! \left(\frac{r}{R}\right)^{2N} U(r)Y_0\,d\tau.$$

**Weak magnetism from electromagnetic data:** For mirror transitions between isospin $T=1/2$ nuclei:

$$b = A\sqrt{\frac{J+1}{J}}\,\mu_{\text{eff}}, \tag{90}$$
where $\mu_{\text{eff}} = \mp(\mu_1 - \mu_2)$ with $\mu_{1,2}$ the magnetic moments of mother and daughter nuclei.

For pure Gamow-Teller transitions within isospin triplets:

$$b^2_\gamma = \eta^6 \frac{\Gamma_{M1}^{iso}}{E_\gamma^3\alpha}, \tag{91}$$

where $\Gamma_{M1}^{iso}$ is the width of the analog isovector M1 gamma transition, $E_\gamma$ is the photon energy, and $\eta = 1$ (or corrected for degeneracy).

### 2. Partially Conserved Axial Current (PCAC)

The Goldberger-Treiman relation relates the pion-nucleon coupling to $g_A$:

$$g_P(q^2) = -g_A(0)\frac{(2M_n)^2}{m_\pi^2 - q^2}. \tag{92}$$

At $q^2=0$, using nucleon and pion masses: **$g_P(0) \approx -229$**. Higher-order corrections (chiral perturbation theory, Gorringe & Fearing 2004; Bhattacharya et al. 2012) modify this at the ~5% level. Experimental measurements from muon capture ($q^2 = 0.88m_\mu^2$) are in sufficient agreement with PCAC results to warrant using them.

**Critical point:** The quenching of $g_P$ due to meson exchange effects in nuclear matter is nucleus-dependent and can reach ~80%. Great care is needed when evaluating pseudoscalar currents for nuclear decays — Eq. (92) applies only to free nucleons.

### 3. Second Class Currents (SCC)

Terms are classified by their transformation under **G-parity**: $G = Ce^{i\pi T_y}$ (charge conjugation × π rotation in isospin space). Two terms transform differently from their vector/axial-vector analogues:
- **Induced scalar** — already eliminated by CVC
- **Induced tensor** — no similar constraint exists

The axial tensor form factor $d$ contains both first-class ($d_I$) and second-class ($d_{II}$) components. Experimentally, **no evidence for SCC has been found** (Minamisono et al., 2001, 2011; Grenacs, 1985), so $d_{II} \equiv 0$. The form factor reduces to $d = d_I$, non-zero only for transitions between **non-analog states**.

---

# Hayen et al. (2017) — Summary of Chunk 8

## Section VI.F: Shake-off and Shake-up Processes; The Endpoint Shift (§6.C, lines 839–1055)

### Overview
Atomic orbital changes due to β decay (screening, exchange) affect the internal "housekeeping" of the atom. Initial and final quantum states are no longer orthogonal, allowing discrete excitations into higher allowed states (**shake-up**) and even continuum ejection (**shake-off**). These reduce available energy for the lepton pair and modify final-state interactions. At the target precision of $10^{-4}$ down to 1 keV, these effects are not always negligible.

---

### 1. Shake-Up (§6.C.1)

#### Mechanism
Atomic electrons in final states can be excited into higher bound states rather than undergoing double ionization (shake-off). The probability is computed from radial wave function overlap integrals of initial and final bound states:
$$P_f = |\langle \phi_n^{0\prime} | \varphi_i \rangle|^2$$
with shake-off probability as the deviation from unity when summing all such probabilities.

**Few-electron systems are particularly susceptible.** In tritium ($^3$H), the single electron in the final state has significant probability to occupy the 2s orbital or higher (Arafune & Watanabe, 1986; Hargrove et al., 1999; Williams & Koonin, 1983).

#### Recoil-Induced Excitation
The recoil momentum of the daughter nucleus induces excitations due to sudden acceleration. The average excitation energy is (Feagin et al., 1979):
$$\Delta E_R = \frac{1}{2} Z m_e v_R^2 = Z E_R \frac{m_e}{M}, \quad \text{(158)}$$

For $Z=50$ and maximum recoil energy of 1 keV, the average excitation energy is ~0.25 eV. Atomic excitations are on the order of tens of eV, so **incomplete overlap between initial/final atomic states dominates** the shake-up picture.

#### Q-Value Reduction
The Q value of the decay is reduced by the mean excitation energy:
$$\Delta E_{ex} = \sum_f P_f E_f, \quad \text{(159)}$$

where $E_f$ is the excitation energy of final state $f$, and $P_f$ is the probability to populate that level.

#### Effect on Screening Correction (Eq. 160)
Possible excitations change screening:
$$S(Z, W) = 1 - \frac{Z\alpha}{p^2} + \sum_h \frac{1}{W} |\langle \phi_n^{0\prime} | \varphi_i \rangle|^2 \left( \sum_{i=1}^n \frac{1}{r_i} \right) = 1 - \frac{Z\alpha}{p^2} + \frac{1}{W} V_{00}, \quad \text{(160)}$$

In tritium β⁻ decay, the screening potential can change by as much as **20%** because the atomic Coulomb interaction involves only a single electron. For systems with more electrons, relative changes decrease but scale with $Z$. Neglecting shake-up introduces errors of order $\sim 1 \times 10^{-4}$.

#### Effect on Exchange Correction
- Medium to high Z: Shake-up probability is negligible at the $10^{-3}$ level.
- Low Z (e.g., tritium): Shake-up probabilities can be significant (~25%), and exchange with higher-lying ns orbitals must be included even if unoccupied in the initial state.

#### Effect on Exchange Correction from Tritium
Inclusion of shake-up changes the exchange correction at the few $10^{-4}$ level down to 0.5 keV for tritium only. Otherwise, negligible.

---

### 2. Shake-Off (§6.C.2)

#### Probability Scaling
Non-orthogonality allows excitations including continuum ejection. Key principles from Carlson et al. (1968):
- **(i)** For given principal quantum number $n$, shake-off probability **decreases** with increasing $Z$ (Pauli principle: only unfilled shells accessible; overlap integrals vanish for inner electrons).
- **(ii)** For a given atom, shake-off probability **increases** with $n$.
- **(iii)** Total shake-off probability is reasonably independent of Z, occurring in ~20–30% of all decays.

K-shell ejection probability: $\sim 2 \times 10^{-3}$ for tritium; a few $10^{-4}$ for higher Z — small effect and not the primary interest here.

#### One-Electron Ionization Probability (Eq. 161)
$$p_\pi^D(M) \simeq 1 - \sum_{n_0 \le n_{\max}} |\langle \phi_{n_0 l'}^{D} | \varphi_{n l}^M \rangle|^2 + K |\langle \phi_{l' \pm 1}^{D} | r | \varphi_{n l}^M \rangle|^2 - K \, \text{Re}[\langle \phi_{l'}^{D} | \varphi_i \rangle \langle \phi_{l'}^{D} | r | \varphi_{n l}^M \rangle], \quad \text{(161)}$$

where $\phi_{nl}^{D(M)}$ has primary quantum numbers $(n,l)$ for daughter/mother, and $K = 2E_R/M$. Holes created by shake-off are filled via X-ray or Auger electron emission.

#### Mean Energy Loss from Shake-Off
For chlorine (Couratin et al., 2013 + Desclaux binding energies), mean energy loss per decay is **~5 eV**. Expected similar magnitude for all Z. To first order, this can be approximated as an additional decrease in $W_0$ by the weighted mean excitation energy.

#### Change in Screening from Shake-Off (Eq. 162–163)
For a hole in outer shell orbital $n,l$:
$$\frac{V_0 - V_{SO}}{V_0} = \frac{\langle f^0 | r^{-1} | i \rangle - \sum f'^0 \langle f^0 | i \rangle \langle f^0 | r | i \rangle}{\langle f^0 | r^{-1} | i \rangle} \approx \frac{1}{\langle f^0 | r^{-1} | i \rangle} \sum_n \frac{\langle f'^0 | i \rangle^2 Z_{\text{eff}}(1 - \Delta Z_{\text{eff}})}{n \cdot n^2}, \quad \text{(162, 163)}$$

With $\Delta Z_{\text{eff}} \approx 0.3$–$0.4$, for Ru$^{1+}$ this yields ~**0.1%** relative change in screening potential → few $10^{-5}$ level in spectral shape (lowest energy range). **Negligible at target precision.**

#### Shake-Off Effect on Exchange
Requires combining Law & Campbell (1972a,b) indistinguishability approach with Harston & Pyper (1992) exchange framework. In medium/high Z nuclei, probability for shake-off to occur in $s_{1/2}$ or $p_{1/2}$ state is **< 0.1%** — effect decreases with increasing Z.

For tritium decay: $\chi_{\text{ex}}^{\text{cont}}$ has maximum value of **-0.01% at 1 keV**. For $^{35}$Ar hole creation in ns states, safely negligible at current precision order.

---

### Section VI.E Continued: Exchange Correction — Analytical Parametrisation (§6.E.4)

#### Challenge
Exchange correction requires radial wave functions of both continuum and bound states over **all space** (not just near origin), for arbitrary potentials. This is a significant challenge for analytical description because:
- Most other effects need wave functions only near the origin (power expansions suffice).
- Exchange needs integration over entire spatial domain.
- Hydrogenic approach with simple screening is insufficient for high precision.
- Pure Coulomb field analytical results [Eq.(156)] do not allow easy insertion of fit parameters.

#### Analytical Fit Formula (Eq. 157)
The exchange correction $X(W)$ as a function of kinetic energy $W_0 \equiv W - 1$:
$$X(W) \approx 1 + \frac{a}{W_0} + \frac{b}{W_0^2} + c \cdot e^{-d W_0} + \frac{e \cdot \sin[(W-f)g+h]}{W^i}, \quad \text{(157)}$$

**9 fit parameters** ($a, b, c, d, e, f, g, h, i$), tabulated in Appendix G for each Z individually. The $Z$ dependence is complex due to atomic shell effects (see Fig. 9).

#### Performance
- Excellent agreement over full tested energy range up to **1 MeV** (where all contributions dip below $10^{-4}$ level).
- Maximum differences: order of a few $\times 10^{-4}$.
- Average residuals: at the **$< 10^{-5}$** level.
- Phase space integral influence: typically **few $\times 10^{-5}$** in absolute terms.

#### Z Dependence — Shell Effects (Fig. 9)
Figure 9 shows orbital exchange contributions for all atoms from $s_{1/2}$ orbitals evaluated at 3 keV. Clear shell closure effects observed:
- **Ar shell closure**: distinct feature around $Z=18$.
- **Kr shell closure**: feature around $Z=36$.
- **Xe shell closure**: feature around $Z=54$.

Increased binding energy at higher Z reduces spatial extension of bound states → decreases overlap between bound and continuum electron radial wave functions at low energies.

---

### Section VI.E Continued: Atomic Potential Sensitivity (§6.E.3, lines 680–791)

Three potentials compared with increasing complexity for exchange correction calculation:

#### 1. Simple Exponentially Screened Field
- Screening strength adjusted for best agreement with bound state energies (decreasing importance for higher $n$).
- Caveat: tuning to match 1s binding energies can result in unbound higher ns states — neglected here since lowest s states dominate.

#### 2. Three Yukawa Potentials
- Coefficients fitted to Dirac-Hartree-Fock-Slater numerical data (Salvat et al., 1987).

#### 3. Complete Potential (Optimized)
- Full potential including exchange optimization turned ON ("Optimized") and OFF ("Unoptimized").
- One free parameter in the exchange potential varied for agreement with numerically calculated binding energies to within **0.1% or 0.1 eV** (Kotochigova et al., 1997a,b; Desclaux, 1973).

#### Results — $^{63}$Ni Comparison (Fig. 8)
- General trend replicated by all potentials apart from very lowest energies where some orbitals give large negative contributions.
- **Without optimization**: binding energies for higher lying orbitals seriously in error → incorrect exchange correction contributions.
- The total net magnitude is a delicate quantity in low energy regions.

#### Nuclear Radius Sensitivity (Fig. 8, bottom)
Varying nuclear radius by 10% with optimized potential:
- Discrepancy grows towards lower energies.
- Only crosses $10^{-4}$ at roughly **~1 keV** for $^{63}$Ni.
- For extreme case of $^{241}$Pu: spectral differences cross $10^{-4}$ level at ~**5 keV**, with constant offset reaching several parts in $10^{-4}$.
- Difference shows reasonably linear behavior on $\Delta R$ — extrapolation to reasonable uncertainties (estimated ~0.03 fm, or ~0.8% for $^{63}$Ni) brings effect below **few $\times 10^{-5}$**.

---

### Numerical Methods: RADIAL Package (§6.E.2, lines 628–677)

#### Wave Function Computation
- Electron radial wave functions $r g_\kappa(r)$, $r f_\kappa(r)$ calculated numerically using **power expansion** at both singular points ($r=0$, $r=\infty$).
- Connected via rescaling of inner solution with $\alpha\kappa$.
- General approach: Behrens & Bühring (1982); Salvat et al. (1995).
- Fortran 77 package **RADIAL** (Salvat et al., 1995) slightly modified and interfaced with custom code for exchange effect calculation.

#### Common Grid Construction
Since free and bound state wave functions are solved on different grids:
- A common grid is constructed.
- Evaluated using exact results or Lagrangian three-point interpolation.
- The overlap integrals $\langle E_{s0} | n s \rangle$ in Eqs.(154)-(155) evaluated directly.
- **Critical**: Grid must be sufficiently dense to avoid systematic errors from interpolation.

#### Example: $^{45}$Ca Exchange Correction (Fig. 6)
- Total effect rises to nearly **10%** in the few-keV region.
- Drops significantly in the 100 eV region due to negative contributions from higher lying orbitals (e.g., 2s).
- At very low energies (~1 keV or lower), contributions from higher-lying orbitals can become **negative**, lowering total correction by several percentage points.

#### Example: $^{45}$Ca p$_{1/2}$ Contribution
- Maximum contribution never exceeds a few parts in $10^5$ — completely negligible for this nucleus.
- Included in all results despite small magnitude.

#### Example: $^{241}$Pu Exchange Correction (Fig. 7)
- Total exchange correction reaches **~29% at 100 eV**.
- p$_{1/2}$ contribution rises to nearly **2% at 100 eV** — definitely not negligible.
- p$_{1/2}$ influence continues over range of tens of keV at target precision.

#### Phase Space Integral Impact
For low endpoint energy transitions, phase space integrals can be significantly altered: up to **30%** for $^{241}$Pu transitions. Effect drops quickly and is typically negligible at several hundreds of keV.

---

### Key Implementation Notes

1. **Exchange correction analytical fit** (Eq. 157) with 9 parameters per Z should be implemented as a lookup table from Appendix G values for fast computation across full energy range.
2. **Atomic potential choice matters critically**: optimized vs unoptimized can differ by $>10^{-4}$ at low energies. The complete potential approach (RADIAL-based) is recommended over simple screened Coulomb models.
3. **Nuclear radius sensitivity**: 10% variation in R causes $>10^{-4}$ spectral differences for heavy nuclei — accurate nuclear charge radii are essential input.
4. **Shake-off/shake-up effects** are generally below target precision except for:
   - Tritium (shake-up ~25%, exchange correction changes at few $10^{-4}$ level).
   - Low-Z light nuclei where shake-up probabilities can be significant.
5. **p$_{1/2}$ orbitals**: While typically ignored, they contribute up to 2% of total exchange in heavy nuclei ($^{241}$Pu) and are non-negligible over tens of keV range — must be included for precision work.

---

# Hayen et al. (2017) — Summary of Chunk 9 (§D Atomic Overlap; §E Bound State β Decay)

## Section VII.D: Atomic Overlap — Alternative Atomic Excitation Correction (lines 1–38)

### Bahcall Correction (Eq. 166)
The β decay results in a sudden change of nuclear potential, so initial and final atomic orbital wave functions only partially overlap. In first approach, this reduces to a difference in atomic binding energies (Bahcall 1963a,c,1965), also included in detailed Ft analyses (Hardy & Towner, 2009).

The correction is constructed by looking at relative change in spectrum shape when changing $W_0$ with $W_0 - \Delta E_{ex}$. Written as:
$$r(Z,W) = 1 - \frac{1}{W_0 - W} \frac{\partial^2 B(G)}{\partial Z^2}, \quad \text{(166)}$$

where $B(G)$ is the total atomic binding energy for a neutral atom with $Z \pm 1$ protons (β∓ decay). The second derivative relates to average excitation energy via $\Delta E_{ex} = -\frac{1}{2}\frac{\partial^2 B(G)}{\partial Z^2}$.

Parametrization using numerical values by Desclaux (1973):
$$K(Z) = -0.872 + 1.270 \, Z^{0.097} + 9.062 \times 10^{-11} \, Z^{4.5}, \quad \text{(169)}$$

The correction $r(Z,W)$ then becomes:
$$r(Z,W) = 1 - \frac{1}{W_0 - W^2} \left[ B(G)'' + 2(C'_0 + C_1) \right], \quad \text{(170)}$$

#### Importance Near Endpoint
- Influence felt mainly near the endpoint where $\frac{\partial^2}{\partial Z^2}B(G)$ can reach **a few hundreds of eV**.
- Particularly important for low energy transitions:
  - **$^{63}$Ni** (endpoint 67.2 keV): correction reaches **1% at 15 keV** before endpoint, increases rapidly past endpoint.
  - **$^{241}$Pu** (endpoint 20.8 keV): similarly important.

---

## Section VII.E: Bound State β Decay (§E, lines 39–67)

### Two-Body Decay Mechanism
Since the electron is created inside an electronic potential well, there exists a possibility for β decay to be captured inside the potential well and produce an electron in a bound atomic state, reducing the decay to a **two-body problem**. First studied by Daudel et al. (1947), expanded by Bahcall (1961), Kabir (1967), Bück (1983), Pyper & Harston (1988).

### Branching Ratio (Eq. 171)
In Bahcall notation, the ratio of bound state to continuum decay probabilities:
$$\frac{\Gamma_b}{\Gamma_c} = \frac{\pi (\alpha Z)^3}{f(Z,W_0)} (W_0 - 1)^2 \Sigma, \quad \text{(171)}$$

where $\Sigma$ depends on binding energy of bound orbital, atomic overlap integrals, and orbital wave functions at nuclear surface.

#### Magnitudes
- **Free neutron decay**: ratio ~ $4.2 \times 10^{-6}$ — completely negligible.
- **Tritium ($^3$H)**: roughly **1% for T⁺** and **0.5% for T** in initial state.
- Higher $W_0$: smaller ratios (phase space integral $f \propto W_0^5$). Kinematic dependence $\Gamma_b/\Gamma_c \propto W_0^{-3}$.
- Superallowed $0^+ \rightarrow 0^+$ decays with Q values of several MeV: probability completely negligible.
- Low energy decays and (partially) ionized initial states: correction can grow significantly.

### Impact on ft Analysis
This does **not** affect the β spectrum shape, as it is a separate final state in S-matrix calculation. It enters the equation when considering the **ft analysis**, however, and cannot always be neglected there.

---

## Section VII.F: Chemical Influences (§F, lines 1–240)

### Overview
In many experiments, the decaying atom is bound within a molecule, modifying electronic structure as electrons rearrange in molecular orbitals. Additional electrons and spectator nuclei influence Coulombic final state interactions; rotational and vibrational states open more possibilities for energy transfer to molecular final state (Cantwell, 1956).

Essential in e.g., antineutrino mass determination in tritium system. Seminal works by Saenz & Froelich describe effects *ab initio* and analytically:
- Doss et al. (2006), Jonsell et al. (1999), Saenz & Froelich (1997a,b), Saenz et al. (2000).

In most precise analysis, this effect introduces a **systematic error of 0.05%**.

---

### 1. Recoil Corrections in Molecules (§F.1)
When β decay occurs inside a molecule:
- The recoiling daughter nucleus moves inside the molecular potential rather than vacuum.
- Potential described by Born-Oppenheimer energy curve, modelled with Lennard-Jones type potential.
- Daughter nucleus kicks molecule into (predissociative) rovibrational state.

**Dissociation possibilities:**
- **Rovibrational states**: total energy > dissociation energy → typical energy scale ~ few eV.
- **Electronic excitation** into resonant continuum: tens to hundreds of eV.

Integrated probability depends on endpoint energy (higher endpoint → higher recoil). At target precision, this effect is neglected as a small correction on already small corrections.

For detecting recoiling daughter atom (e.g., β-ν correlation measurement $a_{\beta\nu}$):
- Probability for dissociation decreases with increasing Z (recoil energy $\propto 1/m$).
- Clear dependence on β-ν angular correlation: higher recoil momenta populate higher rotational/vibrational bands.
- Can be partially included by Monte Carlo simulation of recoiling nucleus in Lennard-Jones potential.

---

### 2. Influence on Q Value (§F.2)
Molecular effects add further decrease to Q value because of **rotational and vibrational excitations**.

Numerical results mainly for $T_2$ (tritium molecule), due to high precision needed for antineutrino mass determination:
- Excitation possibilities change to a broad continuum with resonances.
- Broadening through population of rovibrational states.
- Width of populated rovibrational energy spectrum lies in **few eV region** — typically negligible.
- Excitation into the continuum becomes non-trivial (demonstrated for $T_2$ by Doss et al., 2006).

Averaged over entire continuum: differences between several tritium-substituted molecules is on order of **few eV**, including atomic tritium.

The spectral shape depends on $W_0^2$, so relative error goes as $2 \times \sigma_Q / Q$. For ft values (statistical rate function $\propto Q^5$), precision requirement is even stricter:
- Study case by case for required accuracy.
- Lowest energy transitions (tritium, $^{63}$Ni, $^{241}$Pu): crudely expect **few $10^{-4}$ to $10^{-3}$** level error for $\sigma_Q$ on order of few tens of eV.

Decay rate change treated approximately by Pyper & Harston (1988):
$$\frac{\Delta \lambda}{\lambda} \approx \frac{3 \Delta W_0}{W_0} \left( 1 + \gamma \frac{W_0}{6} \right), \quad \text{(172)}$$

where $\Delta W_0$ is difference in mean endpoint energy after averaging over all final states between two chemical states.

---

### 3. Molecular Screening (§F.3)
Following Saenz & Froelich (1997a): all Coulombic effects treated equally to first order, split into electronic and nuclear parts.

#### Fermi Function First Order (Eq. 173)
$$P_{\text{nuc}}^{(0,1)}(p) = \sum_n \alpha |\langle \phi_0 | \varphi_n \rangle|^2 Z \left[ \frac{8W}{\pi p} + \frac{8p}{3\pi W} + \sum_S Z_S \left( \frac{W}{p^2} + \frac{1}{W} \right) \langle \xi_{000}' | \frac{1}{R_S} | \xi_{000}' \rangle \right],$$

where sum $n$ extends over all final states, sum $S$ takes into account all spectator nuclei, and $R_S$ is distance operator between decaying atom and spectator nucleus S. First part represents Fermi function to first order (larger by factor $p$ than effects of spectator nuclei).

#### Analytical Approximation
For analytical approximation: ignore energy difference in final atomic states; use closure to perform sum over $n$ (first-order corrected using atomic mismatch correction, Sec. VII.C). Each atom considered as having inert inner structure + shared valence wave function.

Full electronic part of molecular wave function (Born-Oppenheimer approximation):
$$\psi_{\text{Mol}} = |\text{valence}\rangle \times \prod_i |\text{inert}_i\rangle, \quad \text{(174)}$$

Using clamped-nuclei approximation for rovibrational ground state: $\langle \xi_{000}' | R_S^{-1} | \xi_{000}' \rangle \approx R_{S,e}$ (equilibrium distance).

#### Total Coulombic Influence (Eq. 175)
$$P_{\text{tot}}^{(0,1)}(p) \approx \alpha \left[ Z\left(\frac{8W}{\pi p} + \frac{8p}{3\pi W}\right) + \frac{W}{p^2} + \frac{1}{W} \sum_S \frac{Z_S - Z_{S,\text{in}}}{R_{S,e}} \right] - \langle i | \frac{1}{r} | i \rangle_0' - Z_{\text{val}} \langle v | \frac{1}{r} | v \rangle,$$

where $Z_S - Z_{S,\text{in}}$ is screened spectator nucleus charge.

#### Molecular Screening Correction (Eq. 176)
$$\Delta S_{\text{Mol}} = \alpha \left[ \frac{W}{p^2} + \frac{1}{W} \sum_S \frac{Z_S - Z_{S,\text{in}}}{R_{S,e}} - Z_{\text{eff}} \langle r^{-1} \rangle_{\text{Val}} \right],$$

where $Z_{\text{eff}} = Z_{\text{Val}} - (Z - Z_{\text{in}})$ and $\langle r^{-1} \rangle_{\text{Val}}$ is average inverse distance of all valence electrons relative to decaying nucleus.

**Example: $^{45}$Ca in CaCl₂ (calcium doubly oxidized, Ca(II))**:
- Bond length measured: 2.437 Å → in natural units $R_{\text{Cl},e} \sim 600$.
- First term $\propto 3 \times 10^{-5}$ — two orders of magnitude smaller than $V_0$.
- Second term: $Z_{\text{eff}} = 2$, effect also on order of **$2 \times 10^{-5}$** with opposite sign.
- Overall molecular deviation from atomic structure: **order of $2 \times 10^{-5}$**.

Since screening correction is typically at (sub)percent level max, molecular deviations expected to have upper limit at **few $10^{-4}$** level.

---

### 4. Molecular Exchange Effect (§F.4)
In a molecule: electronic phase space greatly enlarged and perturbed relative to single atomic state. Electron density near decaying atom can both increase and decrease in molecular bond, so overlap integral changes accordingly.

#### LCAO Approximation
Molecular orbitals constructed from linear combination of atomic orbitals (LCAO). For qualitative argument: approximate that internal orbitals for all atoms are same in molecular as atomic case. Molecular valence orbitals = combination of occupied atomic valence orbitals + energetically close atomic excited states.

- **Ionic bond** (valence electron nearly fully removed): exchange effect $\approx 0$ because internuclear distances >> β Compton wavelength where wave function rapidly oscillates.
- **Reverse extreme** (valence electron density doubled): naively doubles the exchange correction. Conservative treatment: treat last orbital contribution with **100% error bar**.

#### Error Estimate
Taking Fig. 6 as example: introduces **< $1 \times 10^{-4}$** error from 15 keV onwards, grows to **0.5% at 1 keV**. For $^{241}$Pu: valence orbital is 7s (contribution drops below $10^{-4}$ after 3 keV while max quickly grows to 1% in first 0.5 keV).

Effect of wave function change in lattice vs gas studied by Kolos et al. (1988) for tritium — no significant change found. Due to conservative error bar, **completely absorbed and can be neglected**.

---

### Table VII: Complete Overview of β Spectrum Shape Corrections (§VIII start, lines 243–308)

| Item | Effect | Formula | Magnitude |
|------|--------|---------|-----------|
| 1 | Phase space factor | $pW(W_0 - W)^2$ | Unity or larger |
| 2 | Traditional Fermi function | $F_0$ (Eq. 6) | — |
| 3 | Finite size of nucleus | $L_0$ (Eq. 16) | — |
| 4 | Radiative corrections | $R$ (Eq. 47) | — |
| 5 | Shape factor | $C$ (Eqs. 100, 105) | $10^{-1}$–$10^{-2}$ |
| 6 | Atomic exchange | $X$ (Eq. 157) | — |
| 7 | Atomic mismatch | $r$ (Eq. 170) | |
| 8 | Atomic screening | $S$ (Eq. 144)$^a$ | **$10^{-3}$–$10^{-4}$** |
| 9 | Shake-up | See item 7 & Eq. (160)$^b$ | |
| 10 | Shake-off | See items 7, 163, $\chi_{\text{ex}}^{\text{cont}}$(Eq. 164)$^c$ | |
| 11 | Isovector correction | $C_I$ (Eq. 113) | |
| 12 | Distorted Coulomb potential due to recoil | $Q$ (Eq. 45) | **$10^{-3}$–$10^{-4}$** |
| 13 | Diffuse nuclear surface | $U$ (Eqs. 25, 29) | — |
| 14 | Nuclear deformation | DFS (Eq. 40) & $D_C$ (Eq. 135) | — |
| 15 | Recoiling nucleus | $R_N$ (Eq. 41) | — |
| 16 | Molecular screening | $\Delta S_{\text{Mol}}$ (Eq. 176) | — |
| 17 | Molecular exchange | Case by case | — |
| 18 | Bound state β decay | $\Gamma_b/\Gamma_c$ (Eq. 171)$^d$ | Smaller than $1 \times 10^{-4}$ |
| 19 | Neutrino mass | Negligible | < $1 \times 10^{-4}$ |
| 20 | Forbidden decays | Not incorporated | — |

$^a$: Salvat potential of Eq.(147) used with $X$(Eq.145) set to unity.
$^b$: Shake-up effect on screening discussed in Sec.VII.C.1.
$^c$: Case by case evaluation required for shake-off influences on screening/exchange.
$^d$: Does not affect spectral shape but enters ft analysis.

---

### Section VIII: Overview and Crosscheck — Superallowed Fermi Decays (lines 243–270)

#### Comparison with Towner & Hardy (2015) f Values (Fig. 10)
Ratio of $f$ values calculated within this framework vs. results by Hardy & Towner (2015):

- **Excellent agreement**: all residuals in the few $10^{-4}$ region.
- Particularly interesting for heaviest nuclei: e.g., **$^{74}$Rb** — 11 neutrons away from closest stable isotope, subject to strong deformations and shape coexistence ($\beta_2 = +0.401$).
- For four exotic nuclei ($^{62}$Ga, $^{66}$As, $^{70}$Br, $^{74}$Rb): shown results for both spherical shell model and deformed filling (valence nucleons outside $N=Z=28$ shell closure in Nilsson model).
- Deformation does not influence isovector correction as long as $w$ remains equal.

**Conclusion**: Can trust validity of approach for pure Fermi decays in completely analytical manner without additional computation, even for extremely exotic nuclei.

---

### Mirror Decays — ft Value Comparison (lines 271–450)

The ratio of vector and axial vector f values:
$$\equiv \frac{2Ft_{0^+ \rightarrow 0^+}}{Ft_{\text{mirror}}} = \frac{1}{1 + f_{VA} \rho^2}, \quad \text{(177)}$$

where $\rho$ is the mixing ratio. Separate f values calculated by Towner & Hardy (2015).

#### Vector Sector (§VIII, lines 314–332)
- **Exquisite agreement**: all differences smaller than $4 \times 10^{-4}$.
- Holstein form factor $d$ vanishes identically for mirror nuclei.
- CVC allows precise calculation of weak magnetism from experimental data (Severijns et al., TBP).

#### Axial Vector Sector — Two Approaches (§VIII, lines 327–398)

**Approach A: Extreme Single Particle (ESP) with Spherical Harmonic Oscillator (Fig. 11)**
- General agreement good but distinct features for outliers at **$^{33}$Cl and $^{35}$Ar**: disagreement reaches **~1%**.
- Shell model also unable to correctly calculate weak magnetism contribution ($b/Ac_1$) in these cases — accuracy only ~10%.
- Deviations from unity don't necessarily mean failure of $\Lambda$ evaluation.

**Approach B: Deformed Woods-Saxon (DWS) Potential (Fig. 12)**
- Moving from spherical harmonic oscillator to deformed Woods-Saxon while retaining ESP shows **remarkable agreement**.
- Single particle wave function expanded in spherical harmonic oscillator basis; coefficients calculated by numerical diagonalization of deformed Hamiltonian (C++ code, Hayen & Severijns, TBP).
- Deviations when using CVC results can be attributed to shell model limitations (~10% accuracy on $b/Ac_1$).

#### Induced Pseudoscalar Contribution (§VIII, lines 423–434)
- Assuming free nucleon value $g_P = -229$, contribution comparable to ratio of matrix elements → strong deviations expected.
- Possible shifts even on per mille level when ignored.
- Reasonable error must be attributed to all ft$_A$/ft$_V$ calculations until accurate accounting method found.

#### Key Finding: Weak Magnetism Uncertainty (§VIII, lines 402–416)
Differences in weak magnetism predictions can shift $f_A/f_V$ values by **several parts in $10^3$**, undermining claims of reaching 0.01% theoretical precision (Towner & Hardy, 2015).

Most precisely measured mirror ft isotope: **$^{19}$Ne** — behaves well under ESP approximation and not much affected by this concern. Experimental campaigns ongoing for $a_{\beta\nu}$ of $^{32}$Ar, $^{19}$Ne, $^{35}$Ar; $\beta$-asymmetry parameter of $^{35}$Ar, $^{37}$K.

For extraction of $V_{ud}$ to be valid: significant attention needed for precise evaluation of **$\Lambda$ factor in Eq.(107)**.

---

# Hayen et al. (2017) — Summary of Chunks 10, 11, 12 (§VIII–X + Appendices A–G)

---

## Section VIII: Overview and Crosscheck (continued from chunk 9)

### Mirror Decays — Deformed Woods-Saxon Results (Fig. 12)
Moving to deformed Woods-Saxon potential while retaining ESP shows **remarkable agreement** with much more advanced shell model calculations, using only the extreme single particle approximation. The C++ code (Hayen & Severijns, TBP) calculates all nuclear relevant matrix elements $V/A \, M_{KL}^{s}$ automatically.

### Remaining Uncertainties: Weak Magnetism
- Differences in weak magnetism predictions can shift $f_A/f_V$ values by **several parts in $10^3$**.
- Shell model performs only marginally better than deformed ESP approach for predicting weak magnetism (~10% accuracy).
- Raises serious questions on validity and accuracy of presently available calculations at claimed 0.01% precision (Towner & Hardy, 2015).

### Induced Pseudoscalar Component (§VIII end)
Assuming free nucleon value $g_P = -229$, contribution comparable to ratio of matrix elements → strong deviations expected on per mille level when ignored. Reasonable error must be attributed to all ft$_A$/ft$_V$ calculations until accurate accounting method found.

---

## Section IX: Beta-Spectrum Shape Sensitivity to Weak Magnetism and Fierz Terms (§IX, lines 452–71)

### Normalized Shape Factor (Eq. 178)
$$S(W) = C(Z,W) / C(Z, W_{\text{norm}}), \quad \text{(178)}$$

where $C$ is the shape factor normalized to unity at some reference energy $W_{\text{norm}}$. The slope $dS/dW$ provides physics information. Typical value for $b/Ac_1 \approx 5$ (Wauters et al., 2010). For $b_{\text{Fierz}}$: precision of **1% or better** typically required in view of current constraints on CS/CT coupling constants.

At ~0.1% precision, measurements of $b_{\text{Fierz}}$ in nuclear β decay and neutron decay remain competitive with direct searches for new bosons (scalar/tensor type) at LHC ($p + p \rightarrow e^+ + \cancel{E}_T + X$).

### Spectral Sensitivity (§IX.A, lines 1–37 from chunk 11)
- Weak magnetism and non-zero Fierz interference term modify β spectrum shape in energy-dependent fashion.
- **Fierz term**: sensitivity varies as $1/E$ → advantageous for measuring at **low energies**.
- **Weak magnetism**: effect proportional to $E$ → favorable for measuring at **high energy endpoint**.

#### Experimental Shape Factor (historical approach — Calaprice & Holstein, 1976)
Unnormalized shape factor: counts per β-particle energy divided by phase space factor $K(Z,W,W_0,M) \cdot pW(W-W_0)^2$. Then normalized to unity at $W_{\text{norm}}$ by dividing by value at that energy.

#### Slope Analysis
For $b_{\text{Fierz}} = 0.005$, slope over first 250 keV interval is **−0.33% MeV⁻¹**. Energy-dependent information of theoretical spectrum must be accurate enough to guarantee all remaining slope artifacts are **< 0.1% MeV⁻¹**.

#### Low-Energy Degeneracy Problem
At lowest energies, Fierz term spectral modification is approximately linear in energy regime: experiments sensitive to $\sim (b_{\text{wm}} - b_{\text{Fierz}})E$ — **no way to decouple separate contributions** from their energy-dependent behavior.

This underscores the required precision and accuracy in evaluating weak magnetism contribution when searching for new physics.

#### Example: Shape Factor Plot (§IX, Fig. 13)
Shape factor $S(W)$ shown for:
- $b/Ac_1 = 5$ and $b/Ac_1 = 3$ (weak magnetism)
- $b_{\text{Fierz}} = 0.005$ and $b_{\text{Fierz}} = 0.001$

Both effects normalized at 100 keV. Clearly shows need for accurate value of $b/Ac_1$ when extracting Fierz term, and the advantage of measuring low energy transitions for Fierz vs high energy endpoint for weak magnetism.

---

## Section X: Conclusions and Outlook (§X, lines 38–62 from chunk 11)

Key conclusions:
1. **Fully analytical description** of β spectrum shape presented — combining kinematical, electromagnetic, nuclear and atomic corrections to a few parts in $10^4$. Corrections expected to hold from **1 keV to the endpoint**.
2. **Atomic corrections importance**: largest deviations at low energy but influence felt at high energy too. Screening framework combined with Dirac-Hartree-Fock-Slater calculations shows excellent agreement down to lowest energies — avoids ambiguous screening exponent evaluation.
3. Atomic exchange calculations performed for entire atomic number chain; fit provided for analytical evaluation with high precision. **Exchange with p$_{1/2}$ orbitals** (oft-neglected) shown to be significant and included.
4. Combined with Wilkinson's electromagnetic finite size corrections.
5. Nuclear structure corrections derived for both Fermi and Gamow-Teller decay in Holstein notation while maintaining Behrens-Bühring precision of lepton wave function and nuclear current decomposition.
6. **Weak magnetism**: can be analytically calculated to sufficient precision given proper extreme single-particle description. Deformed Woods-Saxon potential provides excellent agreement when spherical HO fails.
7. Induced pseudoscalar contribution influences Gamow-Teller shape factor — important for mirror decay experiments reaching superallowed transition precision.
8. **Compared against Towner & coauthors** most precise numerical f value calculations: very good agreement throughout entire investigated mass range, even for exotic deformed nuclei.
9. Spectral shape measurements have advantage of being less sensitive to nuclear structure effects outside of weak magnetism → basis for experimental BSM searches exploring per mille regime.

### Reactor Antineutrino Anomaly (§X end)
- Beta-spectrum shape calculation is perfect tool to study **weak magnetism contribution**, particularly important for higher masses where little is known experimentally.
- Important ingredient in analysis of reactor antineutrino anomaly (Mueller et al., 2011).
- Current analysis assumes constant weak magnetism throughout entire fission fragment region — not optimal treatment of this complex effect.
- Correct treatment of higher order effects forms basis for correct translation from **electron to antineutrino spectra**.
- Atomic effects, though mainly confined to lower energy regions, are transported to end of antineutrino spectrum and remain highly relevant → underlines importance of precise atomic description at low energy.

### Automation (§X)
Full calculation automated in custom code (Hayen & Severijns, TBP): based on simple configuration files, aims for flexibility and user-friendliness given vast size of work.

---

## Appendix A: General Shape Factor — Charge Distribution Corrections (§A)

### Superallowed Fermi Decay with Full Charge Distribution
Using Wilkinson (1993b) notation to write form factors $F^{(0)}_{ke\,mn\rho}$ and $F^{(1)}_{ke\,mn\rho}$:

$$F^{(0)}_{ke\,mn\rho} = \frac{V F^{(0)}}{V F^{(0)}_{000}}, \quad F^{(1)}_{ke\,mn\rho} = \dots, \quad \text{(A2)}$$

Due to property $I(k_e,m,n,0) \equiv 1$:
$$V F^{(1)}_{000} / V F^{(0)}_{000} = \langle r^2 \rangle / R^2, \quad \text{(A3)}$$

Assuming uniformly charged sphere with radius $R$ → $\langle r^2 \rangle_{\text{exp}} = 3R^2/5$. Better: use **modified Gaussian distribution** (Eq. 18) with one fit parameter $A$, calculated by Wilkinson's method from experimental charge radii.

#### Modified Gaussian Fit Parameters (Eqs. A4–A6)
$$F^{(1)}_{1111} = 0.757 + 0.0069(1 - \exp(-A/1.008)), \quad \text{(A4)}$$
$$F^{(1)}_{1221} = 0.844 - 0.0182(1 - \exp(-A/1.974)), \quad \text{(A5)}$$
$$F^{(1)}_{1222} = 1.219 - 0.0640(1 - \exp(-A/1.550)), \quad \text{(A6)}$$

#### General Form Factors (Eqs. A1a–c for Fermi; A7a–c for Gamow-Teller)
Full expressions given in terms of $V F^{(n)}$ and $A F^{(n)}$ with various $(k,m,n,\rho)$ combinations, including $\alpha ZR$ expansion terms up to order $(\alpha Z)^2$.

### Screening Within Nuclear Volume (§A end)
Neglected atomic electron influence on electrostatic potential in shape factor evaluation. Largest effect already accounted for by screening correction (Eq. 140/144).

For high-Z cases, perturbation order of expressions may be insufficient:
- First-order expansion gives $\Delta U \sim O(\alpha)$.
- Shape factor changes by **few parts in $10^4$** up to large Z — insignificant at that precision level.
- For best performance: further expansion required or replaced with fully numerical approach.

---

## Appendix B: Single-Particle Matrix Elements (§B, lines 207–326)

### Gamow-Teller Form Factors Not Transformable via CVC
In Behrens-Bühring notation: $V F^{(0)}_{111}$ (weak magnetism), $A M^{(0)}_{110}$ (induced tensor), $A M^{(0)}_{121}$.

#### General Expression (Eq. B1)
$$\sqrt{2J_i+1}\, A/V \, M_{1L}^s \approx \sqrt{\frac{2}{2J_i+1}} G_{1L}^s(\kappa_f,\kappa_i) \frac{\langle r^L \rangle}{R} L, \quad \text{(B1)}$$

#### Spin-Reduced Matrix Elements (Eq. B2)
For $A/V \, M_{11s}$: combination of large and small radial functions with $G_{11s}$ coefficients capturing spin-angular information (Weidenmüller 1961; Behrens & Bühring 1982).

$$\frac{A/V}{\sqrt{2J_i+1}} M_{11s}^s \approx \text{sign}(\kappa_i) G_{11s}(\kappa_f,-\kappa_i) \langle g_f | r/R | f_i \rangle + \text{sign}(\kappa_f) G_{11s}(-\kappa_f,\kappa_i) \langle f_f | r/R | g_i \rangle, \quad \text{(B2)}$$

### Explicit Single-Particle Results (categorized by $j_i$, $j_f$)
For $j_i = j_f = l + 1/2$:
$$M_{0101}^0 = \sqrt{\frac{2(l+1)(2l+3)}{(2J_i+1)(2l+1)}} \cdot I, \quad \text{(B3a-type)}$$
$$M_{0121}^0 = -\sqrt{\frac{2}{2J_i+1}} \sqrt{\frac{l(l+1)}{(2l+1)(2l+3)}} \cdot 2\frac{\langle r^2 \rangle}{R^2}, \quad \text{(B3b-type)}$$
$$M_{0111}^0 = -\sqrt{\frac{6(l+1)(2l+3)}{(2J_i+1)(2l+1)}} \cdot \frac{I}{2M_N R}, \quad \text{(B3c-type)}$$
$$M_{0110}^0 = -\sqrt{\frac{2}{2J_i+1}} \sqrt{\frac{2(l+1)}{(2l+1)(2l+3)}} \cdot \xi\frac{\langle r^2 \rangle}{2R}, \quad \text{(B3d-type)}$$

Similar expressions for $j_i = j_f = l - 1/2$ and mixed cases.

#### Small Function Ratio (Eq. B4–B5)
$$I = \int_0^\infty g_f(r) f_i(r) r^2 dr \approx 1 \text{ in non-relativistic approx.}$$
$$\xi \frac{\langle r^2 \rangle}{R^2} = \int_0^\infty g_f(r) \left[ E_i - E_f - (V_i - V_f) \right] \frac{r^2}{R^2} g_i(r) dr, \quad \text{(B4)}$$

For spherical harmonic oscillator:
$$\xi = \frac{2\nu}{M_N} [2(n_i - n_f) + l_i - l_f], \quad \text{(B5)}$$

where $\nu$ is harmonic oscillator parameter (Eq. 115a).

---

## Appendix C: Many-Particle Matrix Elements in jj-Coupling (§C, lines 344–427)

### Odd-Z, Odd-N to Even-Z, Even-N β Decays and Vice Versa
Using Rose & Osborn (1954a), written using isospin formalism (Wilkinson, 1969). For two particles in initial/final states coupled to core isospin.

#### o-o to e-e Transitions (Eq. C1)
$$\langle j_2 j_2 J_f M_f T_f T_{3f} | \sum_n O_{KLs}^n \tau_n^\pm | j_1 j_2 J_i M_i T_i T_{3i} \rangle = \hat{J}_i \hat{J}_f \hat{T}_i \hat{T}_f (-)^{\dots} \begin{Bmatrix} T_f & 1/2 \\ T_i & T_f^{(2)} \end{Bmatrix} \times [\dots] \langle j_2 || O_{KLs} || j_1 \rangle + [1 - (-)^{j_1+j_2}] \begin{Bmatrix} j_2 & J_f & j_1 \\ J_i & j_2 & K \end{Bmatrix} \langle j_2 || O_{KLs} || j_2 \rangle \delta_{j_1 j_2} \langle 1/2 || t || 1/2 \rangle$$

For e-e to o-o: equivalent formula. $\langle 1/2 || t || 1/2 \rangle = \sqrt{3/2}$.

### Deformation Case
When deformation present, angular momentum $J$ is no longer good quantum number; Eq.(C1) must be rewritten using spin-reduced matrix elements (Berthier & Lipnik, 1966):

$$\langle \phi(J_f K_f; \mathcal{F}_f) || \sum_n O_{KLs}^\pm \tau_n^\pm || \phi(J_i K_i=0; \mathcal{F}_i=0) \rangle = (-)^{\dots} [1 + (-)^{J_i}] \sum_{j_2,j_1} C_{j_2\mathcal{F}_2} C_{j_1\mathcal{F}_1} (-)^{j_2-\mathcal{F}_2} \langle j_2 || O_{KLs} || j_1 \rangle, \quad \text{(C2)}$$

---

## Appendix D: Relativistic Coulomb Amplitudes (§D)

### Smallness of Relativistic Matrix Elements
For spectral shape measurements at aimed-for precision level, relativistic matrix elements are **insignificant** but discussed briefly for completeness.

#### Pure Fermi Transitions (Eq. D1)
$$V C(Z,W)_{\text{rel}} = V F^{(0)}_{011} f_2(W) + V F^{(0)}_{011} f_3(W), \quad \text{(D1)}$$

Coulomb functions $f_2$ and $f_3$ written in terms of Table VIII elements:
$$f_2^V(W) = -2(D_1 + N_1) + 2\frac{\mu}{W} d_1, \quad \text{(D2a)}$$
$$f_3^V(W) = -2(D_3 + N_1 H_2 - N_2 D_1 - N_3) + 2\frac{\mu}{W}(d_3 - N_2 d_1), \quad \text{(D2b)}$$

Note: $D_1(1) + N_1(1) = \alpha Z/2 + W_0 R/3$ — independent of $W$. Typical values for $f_i(W)$ are on percent level smaller. Relativistic matrix elements $\propto V F^{(N)}_{011}$ obey recursion relations after invoking CVC (Eq. 89).

### Table VIII: Coefficients in Leptonic Radial Wave Function Expansion
Coefficients $H_{2\sigma}$, $D_{2\sigma+1}$, $h_{2\sigma}$, $N_\sigma$, and $d_{2\sigma+1}$ for expansion of leptonic radial wave functions using older Behrens-Jänecke (1969) leptonic wave expansion.

---

## Appendix E: Comparison — Behrens-Bühring vs. Holstein Formalisms (§E, lines 462–738)

### Motivation
Confusion exists in the name "finite size effects" as it entails different things for different authors. This work based on rigorous Behrens-Bühring approach but written in more transparent Holstein form factors for clarity and interpretation of experimental results.

---

### 1. Generalization Including Electrostatics (§E.1)
#### Simple Hamiltonian (Eq. E1)
$$H_\beta = \frac{G \cos\theta_C}{\sqrt{2}} [\bar{u}(p)\gamma_\mu(1+\gamma_5)v(l)] \langle p_2 | V^\mu + A^\mu | 0 \rangle + H.c., \quad \text{(E1)}$$

Plane waves used for leptons. S-matrix developed to first order. Coulomb interaction must be incorporated between nucleus and outgoing leptons.

#### Two Approaches:
**Behrens-Bühring (BB):** Expands electron wave function in powers of $(\alpha Z)\rho$, $(WR)^{\nu-\rho}$, $(m_e R)^{\mu-\nu}$ with coefficients encoded in functions $I(|\kappa|,\mu,\nu,\rho; r)$ — sensitive to nuclear potential.

**Holstein et al.:** Starts from generalized matrix element (Armstrong & Kim 1972b):
$$M = \sqrt{\frac{G \cos\theta_C}{2}} \int d^3r\, \bar{\Psi}_e(r,p) \gamma_\mu(1+\gamma_5) v(l) \int \frac{d^3k}{(2\pi)^3} e^{i\mathbf{r}\cdot\mathbf{k}} \langle p_2 + p - k | V^\mu + A^\mu | 0 \rangle, \quad \text{(E2)}$$

Both results identical up to this point after neglecting difference in initial/final Coulomb potentials. Both expand nuclear current using form factors (replacing impulse-approximation results of Eqs. 83).

---

### 2. Approximations by Holstein and Notes for the Wary (§E.2)
#### Electron Wave Function Expansion Difference
Holstein expansion keeps only leading $j=1/2$ and $j=3/2$ terms → function of $f_\kappa, g_\kappa$ with $\kappa \in \{-2,-1,1,2\}$.

**First critical difference**: Holstein expansion performed **only for uniformly charged sphere**. Discussion in Sec. IV shows this is insufficient — diffuse charge distribution effects cannot be ignored.

#### Nuclear Current Expansion Difference
Second crucial difference: combination of nuclear current expansion with lepton wave function expansion handled differently.

---

### 3. Fermi Function and Finite Size Correction (§E.3)
Coulomb corrections depend on integrals A–G, which depend on weak charge distribution (Calaprice & Holstein, 1976).

#### Gamow-Teller Approximation: Constant Form Factor (§E.3)
If $c(q^2)$ approximated as constant → $\rho(r)$ becomes Dirac delta function → only A and B survive:
$$|A|^2 + |B|^2 = N_e \frac{1}{2} \frac{2m_e}{W} \int d^3r\, \delta^3(r) [g_{-1}^2(r) + (W_0 - W)^2/p^2 f_1^2(r)] = F_0 L_0, \quad \text{(E9)}$$

Taking into account differing normalization definitions of BB vs Holstein.

Holstein result by Calaprice & Holstein (1976) and Huffaker & Laird (1967):
$$|A|^2 + |B|^2 \approx F_0 [1 \mp \alpha Z W R / 15], \quad \text{(E10)}$$

**Clearly less precise than Eq.(16)** — contains no corrections from diffuse nuclear charge, for which explicit expressions and analytical parametrizations derived in Sec. IV.A.2.

#### Nuclear-Leptonic Convolution (§E.4)
Finite nuclear size affects not only Coulomb potential but behavior of wave function near origin. Matrix element requires averaging over all nucleon positions (inherent in C factor definition, Sec. VI.E).

First-order expansion: $c(q^2) \approx c_1 + c_2 q^2$:
$$\rho \approx \int \frac{d^3k}{(2\pi)^3} e^{i\mathbf{r}\cdot\mathbf{k}} [1 + (W_0 + i\nabla)^2] = [1 + (W_0 + i\nabla)^2] \delta^3(r), \quad \text{(E12)}$$

All leptonic radial wave functions evaluated at nuclear center, multiplied by constant factor depending on nuclear structure → all results proportional to Fermi function defined at origin.

#### Critical Comparison: Holstein vs BB Results (§E.4)
For uniformly charged sphere with $c_2/c_1 = \frac{1}{6} R^2$:
$$HS\, C(Z,W)^0_{\text{old}} = 1 + \frac{R^2}{5} - \frac{(W_0 R)^2}{5} - \frac{(\alpha Z)^2}{20} \pm \frac{\alpha Z W_0 R}{15} \mp \frac{\alpha Z W R}{15} + \dots, \quad \text{(E14)}$$

With old expansion (Behrens & Jänecke 1969; Wilkinson 1990) for Gamow-Teller:
$$old\, C(Z,W)^0 = 1 + {}^A C_0 + A C_1 W + A C_2 W^2, \quad \text{(E15)}$$

where ${}^A C_0$, $C_1$, $C_2$ given by Eq.(E16) — every term identifiable with similar one in Eq.(E14), **except factor 2 difference** in $\alpha Z W_0 R$ term. All results valid only for uniformly charged sphere.

In Sec.A: deviation from uniform density analyzed; isovector correction $C_I$ (Sec.VI.F) introduced to correct breakdown — due to its large contribution, cannot be neglected.

---

### 4. Induced Terms (§E.5)
Holstein formalism requires adjustment for approximation in Eq.(E5). Initially, Coulomb corrections to induced terms neglected (already suppressed by $q/M$), but not strictly valid per Bottino et al. (1973, 1974).

Starting from Armstrong & Kim result: $\langle f(p_f) | \rangle$ should be replaced by $\langle f(p_f + p_e - p') |$, introducing Coulomb corrections to induced terms as function of $q' = (p_f + p_e - p') - p_i$ rather than $q = p_f - p_i$.

The correction factor from Calaprice & Holstein:
$$\delta h_1(Z) \approx \sqrt{\frac{10}{6}} \frac{\alpha Z}{MR} c_1 (2b + d_{II} \pm d_I \pm c_1), \quad \text{(E18)}$$

Can be obtained from Bottino et al. results after proper conversion. **No artificial separation** of nuclear structure and Coulombic terms in BB formalism → these corrections occur naturally in C factor.

---

### 5. Validity of Harmonic Oscillator Functions (§F)
Wang et al. (2016) addressed finite size corrections using Density Functional Theory (DFT), calculating weak charge density analytically to first order in $\alpha Z$ and R, rewritten in terms of Zemach moments:

$$\delta_{FS}^{\text{Wang}} = -\frac{\alpha Z}{3} \left[ 4 W \langle r \rangle^{(2)} + W_r \langle r \rangle^{(2)}_r - \frac{1}{W} (2\langle r \rangle^{(2)} - \langle r \rangle_r^{(2)}) \right], \quad \text{(F1)}$$

where $\langle r \rangle^{(2)}$ and $\langle r \rangle_r^{(2)}$ defined by Eqs.(F2)–(F3):
$$\langle r \rangle^{(2)} = \int d^3s \int d^3r\, \rho_w(r)\, \rho_{ch}(|\mathbf{r}-\mathbf{s}|), \quad \text{(F2)}$$
$$\langle r \rangle_r^{(2)} = \int d^3s \int d^3r\, \rho_w(r)\, s \frac{\partial}{\partial r} \rho_{ch}(|\mathbf{r}-\mathbf{s}|). \quad \text{(F3)}$$

Here $\rho_w$ and $\rho_{ch}$ are weak and regular charge densities.

#### Comparison: DFT vs Analytical Harmonic Oscillator (§F end)
Table IX compares $hr^2 i_w^{\text{Wang}}$, $hr^2 i_w^{\text{HO}}$, $hri^{(2)}_{\text{BB}}^U$, $hri^{(2)}_{\text{BB}}^{HO}$, and $hri^{(2)}_{\text{Wang}}$ for several β⁻ transitions ($A=14$ to $A=121$).

**Approximations used (Eqs. F4–F5):**
$$\langle r \rangle^{(2)} \approx \langle r \rangle_{ch} + \frac{\langle r^2 \rangle_w \langle r^{-1} \rangle_{ch}}{3}, \quad \text{(F4)}$$
$$\langle r \rangle_r^{(2)} \approx \frac{2\langle r^2 \rangle_w \langle r^{-1} \rangle_{ch}}{3}. \quad \text{(F5)}$$

Eq.(F5) is a **poor approximation** — overestimates integral by ~30% (see Table IX last column for Sn/Sb).

#### Agreement Assessment
- Lower masses: fair agreement.
- Highest masses ($A > 100$): differences appear due to two reasons:
  1. Considering only final decaying neutron not valid — quantum number difference between proton/neutron orbitals too large; actual overlap integral requires $\int r^2 R_{nl}^\alpha R_{n'l'}^\beta dr$.
  2. HO functions are charge-insensitive, while DFT results are not → expect additional effect from isospin invariance breakdown (nuclear mismatch $\delta_C$, Towner & Hardy, 2010).

Application #1 yields results accurate to **within 10% of DFT results**. Remaining discrepancy attributed to inherent limitations of HO treatment and nuclear mismatch problem. The latter drops out to first order when taking ratios of matrix elements suffering from same issue — sufficient for current purposes. Extension to deformed nuclei straightforward (deformed Woods-Saxon with HO basis functions → excellent agreement).

---

## Appendix G: Tabulation of Fit Parameters for Exchange Correction (§G)

### Table X: Eq.(157) Coefficients for Z = 2 to Z = 120
Full tabulation provided in original paper (lines 834–957 of chunk 12). The analytical fit parameters $a, b, c, d, e, f, g, h, i$ from:

$$X(W) \approx 1 + \frac{a}{W_0} + \frac{b}{W_0^2} + c \cdot e^{-d W_0} + \frac{e \cdot \sin[(W-f)g+h]}{W^i}, \quad \text{(157)}$$

where $W_0 = W - 1$ is the kinetic energy. These parameters are already present in the codebase at `/home/unzhakov/workspace/beta_shape_pnpi/beta_spectrum/data/exchange_coeff.csv`.

The table shows clear shell-structure effects (e.g., parameter discontinuities near Ar, Kr, Xe closures) due to atomic orbital occupancy changes.

---

## Summary of All Chunks Read

| Chunk | Content | File |
|-------|---------|------|
| 1–6 | Abstract through nuclear structure intro; all major correction terms | `/tmp/hayen2017_summary_01-06.md` |
| 7 | β decay Hamiltonian, BB vs HS notation, CVC/PCAC/SCC | `/tmp/hayen2017_summary_07.md` |
| 8 | Isospin breakdown; exchange correction (analytical fit §E.4, potential sensitivity §E.3, RADIAL package); shake-off/shake-up (§C) | `/tmp/hayen2017_summary_08.md` |
| 9 | Atomic overlap/Bahcall correction; bound state β decay; chemical influences (§D–F); Table VII overview; superallowed Fermi crosscheck; mirror decays | `/tmp/hayen2017_summary_09.md` |
| 10 | Mirror decays (cont.); Fierz/weak magnetism sensitivity intro; Section VIII end with Table VII complete; BSM physics outlook (§IX start) | part of chunk 10 in summary above |
| 11 | Beta-spectrum shape sensitivity to weak magnetism and Fierz terms (§IX); Conclusions & Outlook (§X) | `/tmp/hayen2017_summary_11.md` (to be written — content covered above) |
| 12 | Appendices A–G: general shape factor; single-particle matrix elements; many-particle jj-coupling; relativistic Coulomb amplitudes; BB vs HS comparison; HO validity; exchange fit parameter table X | part of chunk 12 in summary above |
