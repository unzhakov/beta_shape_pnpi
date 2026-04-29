# Algorithm for Radiative Corrections to Beta-Decay Spectra

## 1. Physical Overview

The radiatively corrected beta-decay spectrum takes the form:

\[\frac{d\Gamma}{dW} = \frac{d\Gamma_0}{dW} \cdot R(W, W_0)\]

where \(d\Gamma_0/dW\) is the uncorrected spectrum and \(R(W, W_0)\) is a multiplicative correction factor. The correction is organized by powers of the fine-structure constant \(\alpha\) and the nuclear charge \(Z\):

\[R(W, W_0) = 1 + \delta^{(1)}(W, W_0) + \delta^{(2)}(Z, W) + \delta^{(3)}(Z, W) + \cdots\]

**Important distinction:** Inner radiative corrections (\(\Delta_R^V\), \(\Delta_R^A\)) are energy-independent constants absorbed into effective coupling constants. They are **not** part of the spectrum shape calculation. Only the outer (energy-dependent) correction \(\delta_R(W, W_0)\) modulates the spectral shape.

---

## 2. Input Parameters

| Symbol | Description | Units |
|--------|-------------|-------|
| \(W\) | Total electron energy | \(m_e c^2\) |
| \(W_0\) | Endpoint energy | \(m_e c^2\) |
| \(Z\) | Nuclear charge of daughter | dimensionless |
| \(m_p\) | Proton mass | \(m_e c^2\) (\(\approx 1836.15\)) |
| \(\alpha\) | Fine-structure constant | \(\approx 1/137.036\) |

Derived quantities:
- Momentum: \(p = \sqrt{W^2 - 1}\)
- Velocity: \(\beta = p / W\)
- Antineutrino energy (for neutrino spectrum corrections): \(\hat{W} = W_0 - W\)

---

## 3. Order \(\alpha\): The Universal Sirlin Function \(g(W, W_0)\)

### 3.1 Mathematical Form

The dominant correction, derived in Sirlin (1967) and refined in Sirlin & Zucchini (1986) and Sirlin (1987), is universal—identical for Fermi and Gamow-Teller transitions, electrons and positrons:

\[g(W, W_0) = 3\ln\left(\frac{m_p}{m_e}\right) - \frac{3}{4} + 4\left(\frac{\tanh^{-1}\beta}{\beta} - 1\right)\left[\frac{W_0 - W}{3W} - \frac{3}{2} + \ln\left(\frac{2(W_0 - W)}{m_e}\right)\right] + \frac{4}{\beta}L\left(\frac{2\beta}{1+\beta}\right) + \frac{\tanh^{-1}\beta}{\beta}\left[2(1+\beta^2) + \frac{(W_0-W)^2}{6W^2} - 4\tanh^{-1}\beta\right]\]

The correction to the spectrum is:

\[\delta^{(1)}(W, W_0) = \frac{\alpha}{2\pi} \, g(W, W_0)\]

### 3.2 The Dilogarithm (Spence Function)

The function \(L(x)\) appearing above is the Spence function (dilogarithm):

\[L(x) = \int_0^x \frac{\ln(1-t)}{t} dt = -\text{Li}_2(-x/(1-x)) \quad \text{or equivalently} \quad L(x) = -\text{Li}_2(x)\]

Using the standard dilogarithm \(\text{Li}_2(x) = -\int_0^x \frac{\ln(1-t)}{t} dt = \sum_{k=1}^\infty \frac{x^k}{k^2}\):

\[L\left(\frac{2\beta}{1+\beta}\right) = -\text{Li}_2\left(\frac{2\beta}{1+\beta}\right)\]

### 3.3 Numerical Stability Considerations

**Inverse hyperbolic tangent:**

\[\tanh^{-1}\beta = \frac{1}{2}\ln\left(\frac{1+\beta}{1-\beta}\right)\]

When \(\beta \to 1\) (ultra-relativistic limit), compute using the identity:

\[\tanh^{-1}\beta = \ln\left(\frac{W+p}{m_e}\right) = \ln(W + p) \quad (\text{in } m_e \text{ units where } m_e = 1)\]

For \(\beta \to 0\) (non-relativistic limit), use the Taylor expansion:

\[\frac{\tanh^{-1}\beta}{\beta} = 1 + \frac{\beta^2}{3} + \frac{\beta^4}{5} + \cdots\]

Avoid computing \(\tanh^{-1}\beta\) directly when \(\beta\) is very small.

**Terms containing \(1/\beta\):** When \(\beta < \epsilon\) (e.g., \(\epsilon = 10^{-6}\)), evaluate the expression using a small-\(\beta\) series expansion to avoid division by zero.

**Dilogarithm evaluation:** Use a robust numerical implementation. For \(x > 1\), apply the identity:

\[\text{Li}_2(x) = \frac{\pi^2}{3} - \frac{1}{2}\ln^2(x) - \text{Li}_2(1/x) - i\pi\ln(x) \quad (x > 1)\]

The imaginary part cancels in the real-valued total correction.

### 3.4 Endpoint Behavior and Soft-Photon Resummation

**Critical observation:** The function \(g(W, W_0)\) contains a term proportional to \(\ln(W_0 - W)\), which diverges as \(W \to W_0\). This divergence is physical—arising from the emission of arbitrarily soft real photons—and the total decay rate integrated over all photon energies remains finite.

For the differential spectrum shape, the leading logarithmic behavior near the endpoint is:

\[g(W, W_0) \sim 4\left(\frac{\tanh^{-1}\beta}{\beta} - 1\right) \ln(W_0 - W) + \text{regular terms}\]

**Resummation prescription (Sirlin 1987, 2011):** Replace the logarithmically divergent term with an exponentiated form:

\[\ln(W_0 - W) \cdot t(\beta) \quad \longrightarrow \quad (W_0 - W)^{t(\beta)} - 1\]

where the exponent is:

\[t(\beta) = \frac{2\alpha}{\pi}\left[\frac{\tanh^{-1}\beta}{\beta} - 1\right]\]

**Implementation strategy:**

1. Compute the full \(g(W, W_0)\) including the \(\ln(W_0 - W)\) term.
2. For \(W_0 - W < \Delta_{\text{cut}}\) (suggested \(\Delta_{\text{cut}} \sim 10^{-3}\) in units of \(m_e c^2\)), extract the coefficient of \(\ln(W_0 - W)\):
   
   \[C(\beta) = 4\left(\frac{\tanh^{-1}\beta}{\beta} - 1\right)\]
   
3. Replace \(C(\beta) \cdot \ln(W_0 - W)\) with \(C(\beta) \cdot \frac{(W_0-W)^{t(\beta)} - 1}{t(\beta)}\) using:
   
   \[\lim_{t \to 0} \frac{x^t - 1}{t} = \ln x\]

4. For energies well below the endpoint (\(W_0 - W \gg \Delta_{\text{cut}}\)), no resummation is needed and the original logarithmic form suffices.

**For low endpoint energies** (e.g., tritium with \(W_0 \approx 1.02\)), the endpoint resummation has negligible effect and can be omitted for simplicity.

---

## 4. Order \(Z\alpha^2\): Finite Nuclear Size Corrections

### 4.1 Structure

The \(O(Z\alpha^2)\) correction, derived in Sirlin (1987), depends on the nuclear charge distribution:

\[\delta^{(2)}(Z, W) = Z\alpha^2 \sum_{i=1}^4 \Delta_i(W)\]

The four contributions \(\Delta_i\) arise from distinct classes of Feynman diagrams (see Sirlin 1987, Sec. 2). The dominant energy dependence comes from \(\Delta_1^0(W) + \Delta_4(W)\), for which the relativistic limit gives:

\[\Delta_1^0(W) + \Delta_4(W) = \ln\left(\frac{M}{m_e}\right) - \frac{5}{12}\ln\left(\frac{2W}{m_e}\right) + \frac{43}{18}\]

where \(M\) is the nucleon mass (use proton mass \(m_p\)).

### 4.2 Nuclear-Model-Dependent Part

The finite-nuclear-size contribution \(\Delta^F = \Delta_1^F + \Delta_2 + \Delta_3\) depends on the nuclear charge distribution. For the **modified Gaussian model** (Sirlin 1987, Eq. 10):

\[\Delta^F(Z) = \ln\left(\frac{\Lambda}{M}\right) - \kappa_1(Z) + O(\Lambda/M)\]

where:
- \(\Lambda = \sqrt{6}/a\), with \(a\) the rms charge radius
- \(M/\Lambda = r_0 A^{1/3} / 0.665\) (using Wilkinson's \(r_0\))
- \(\kappa_1(Z) = \frac{1}{2}\left[\gamma + \ln(3/2k^2) + 2\alpha/(2+3\alpha)\right]\), with \(\gamma = 0.5772\), \(\alpha = (Z-2)/3\), and \(k^2 = \frac{3}{2}\frac{2+5\alpha}{2+3\alpha}\)

### 4.3 Practical Simplification

For light and medium nuclei (\(Z \lesssim 30\)), the total \(O(Z\alpha^2)\) correction can be tabulated once per isotope as an energy-independent shift (the energy dependence from \(\Delta_1^0 + \Delta_4\) is small compared to the order \(\alpha\) term). **For spectrum shape calculations:**

1. Calculate the energy-dependent part \(\Delta_1^0(W) + \Delta_4(W)\) using the relativistic formula above.
2. Add the energy-independent nuclear-structure part \(\Delta^F(Z)\) from precomputed values (see Sirlin 1987, Table II).
3. If high precision is required for low-energy electrons \((p \ll m_e)\), use the full non-relativistic expressions from Sirlin 1987.

For a **uniformly charged sphere**, \(\kappa_2 = \gamma - 4/3 + \ln(\sqrt{10}) = 0.395\), giving results within ~2% of the modified Gaussian model.

---

## 5. Order \(Z^2\alpha^3\): Higher-Order Corrections

The heuristic estimate from Sirlin (1986, 1987) is:

\[\delta^{(3)}(Z, W) \approx Z^2\alpha^3 \left[a\ln\left(\frac{\Lambda}{W}\right) + b f(W) + \frac{4\pi}{3}g(W) - 0.649\ln(2W_0)\right]\]

with:
\[a = \frac{\pi}{3} - \frac{\sqrt{3}}{2}, \quad b = \frac{4}{3\pi^4}(\pi^2 - \gamma_E) - \frac{\pi^2}{18}\]

And:
\[f(W) = \ln(2W_0/m_e) - 5/6\]
\[g(W) = \frac{1}{2}\left[\ln^2(R m_e) - \ln^2(2W/m_e)\right] + \frac{5}{4}\ln(2RW)\]

where \(R = \sqrt{5/4}\,a\) and \(a\) is the rms charge radius.

**Practical note:** For \(Z \lesssim 30\), \(\delta^{(3)} < 0.01\%\) and can typically be neglected for spectrum shape. For superallowed Fermi transitions with heavy nuclei, this correction contributes ~0.05–0.12% to the total correction (Sirlin 1987, Table III).

---

## 6. Neutrino Spectrum Correction

For the antineutrino/neutrino spectrum (Sirlin 2011, Eq. 11), define:

\[\hat{W} \equiv W_0 - W_\nu, \quad \hat{p} \equiv \sqrt{\hat{W}^2 - 1}, \quad \hat{\beta} \equiv \hat{p}/\hat{W}\]

The correction function \(h(\hat{W}, W_0)\) is:

\[h(\hat{W}, W_0) = 3\ln\left(\frac{m_p}{m_e}\right) + \frac{23}{4} - \frac{8}{\hat{\beta}}\text{Li}_2\left(\frac{2\hat{\beta}}{1+\hat{\beta}}\right) + 8\left(\frac{\tanh^{-1}\hat{\beta}}{\hat{\beta}} - 1\right)\ln\left(\frac{2\hat{W}\hat{\beta}}{m_e}\right) + 4\frac{\tanh^{-1}\hat{\beta}}{\hat{\beta}}\left[\frac{7+3\hat{\beta}^2}{8} - 2\tanh^{-1}\hat{\beta}\right]\]

The neutrino spectrum correction is:

\[\delta_\nu(\hat{W}, W_0) = \frac{\alpha}{2\pi} \, h(\hat{W}, W_0)\]

**Key difference from the electron correction:** The neutrino correction \(h\) is finite in the \(m_e \to 0\) limit, while the electron correction \(g\) diverges logarithmically. The \(m_e \to 0\) limit of \(h\) is particularly simple:

\[h(\hat{W}, W_0) \xrightarrow{m_e \to 0} 3\ln\left(\frac{m_p}{2W_0}\right) + \frac{23}{4} - \frac{4\pi^2}{3} - 3\ln\left(\frac{W_\nu}{W_0}\right)\]

This finite limit is explained theoretically by the Kinoshita-Lee-Nauenberg theorem: summing over the degenerate collinear \(e^-\gamma\) final states (as is done when integrating over electron momenta to obtain the neutrino spectrum) cancels the mass singularities (Sirlin 2011, Sec. II.F).

---

## 7. Conversion Between Electron and Neutrino Spectra

Given the electron spectrum \(f_e(W, W_0)\) for a specific decay, the corresponding neutrino spectrum including radiative corrections is (Sirlin 2011, Eq. 19):

\[f_\nu(W_\nu, W_0) = f_e(\hat{W}, W_0) \left[1 + \frac{\alpha}{2\pi}\left(h(\hat{W}, W_0) - g(\hat{W}, W_0)\right)\right]\]

where \(\hat{W} = W_0 - W_\nu\).

The difference function is given explicitly in Sirlin 2011, Eq. 20:

\[h(\hat{W}, W_0) - g(\hat{W}, W_0) = \frac{13}{2} - \frac{4}{\hat{\beta}}\text{Li}_2\left(\frac{2\hat{\beta}}{1+\hat{\beta}}\right) + 4\left(\frac{\tanh^{-1}\hat{\beta}}{\hat{\beta}} - 1\right) \times \Bigg[\frac{11}{6} + \ln\left(\frac{2W_0}{m_e}\right) + 2\ln\hat{\beta} + 2\ln\left(\frac{W_\nu}{W_0}\right) - \ln\left(\frac{W_\nu}{W_0}\right) - \frac{1}{3W_\nu/W_0}\Bigg] + \frac{\tanh^{-1}\hat{\beta}}{\hat{\beta}}\left[\frac{3-\hat{\beta}^2}{2} - 4\tanh^{-1}\hat{\beta} - \frac{(W_\nu/W_0)^2}{6(1-W_\nu/W_0)^2}\right]\]

This conversion is important for reactor neutrino oscillation studies.

---

## 8. Complete Calculation Algorithm

### 8.1 Full Correction Factor for Electron Spectrum

```
Algorithm: RadiativeCorrection(W, W0, Z)
Input:  W  — array of electron total energies [m_e c^2]
        W0 — endpoint energy [m_e c^2]
        Z  — nuclear charge of daughter nucleus
Output: R  — multiplicative correction factor at each W

Constants:
    alpha     = 1/137.035999084
    m_p       = 938.272 / 0.51099895   # proton mass in m_e units
    gamma_E   = 0.5772156649
    pi        = 3.141592653589793

For each W_i in W:
    1. Compute derived quantities:
       p   = sqrt(W_i^2 - 1)
       beta = p / W_i
       
    2. Order α correction (Section 3):
       a. Calculate tanh^{-1}(beta) / beta
          If beta < 1e-6: use Taylor series 1 + beta^2/3 + beta^4/5
          Else: tanh_inv_beta_beta = 0.5 * ln((1+beta)/(1-beta)) / beta
       
       b. Calculate dilogarithm L(2*beta/(1+beta)) = -Li_2(2*beta/(1+beta))
          Use robust numerical implementation
       
       c. Compute g(W_i, W0) per Eq. (20b) of Sirlin 1967
       
       d. If use_endpoint_resummation and (W0 - W_i) < delta_cut:
          t_beta = (2*alpha/pi) * (tanh_inv_beta_beta - 1)
          C_beta = 4 * (tanh_inv_beta_beta - 1)
          Replace: C_beta * ln(W0 - W_i) 
                -> C_beta * ((W0 - W_i)^t_beta - 1) / t_beta
       
       e. delta_1 = (alpha / (2*pi)) * g(W_i, W0)
       
    3. Order Zα² correction (Section 4):
       a. Energy-dependent part:
          Delta_1_4 = ln(m_p) - (5/12)*ln(2*W_i) + 43/18
       
       b. Energy-independent nuclear part:
          Look up Delta_F for the given (Z, A) from precomputed table
          (Use Sirlin 1987 Table II or compute from nuclear model)
       
       c. delta_2 = Z * alpha^2 * (Delta_1_4 + Delta_F)
       
    4. Order Z²α³ correction (Section 5, optional):
       For Z >= 20 or precision better than 0.01% required:
          Compute delta_3 from Sirlin heuristic formula
       Else:
          delta_3 = 0
       
    5. Combine corrections:
       R_i = 1 + delta_1 + delta_2 + delta_3
    
    Return R array
```

### 8.2 Neutrino Spectrum Correction

```
Algorithm: NeutrinoRadiativeCorrection(W_nu, W0, Z)
Input:  W_nu — array of neutrino energies [m_e c^2]
        W0   — endpoint energy [m_e c^2]
        Z    — nuclear charge
Output: R_nu — multiplicative correction factor

For each W_nu in array:
    W_hat = W0 - W_nu
    p_hat = sqrt(W_hat^2 - 1)
    beta_hat = p_hat / W_hat
    
    Compute h(W_hat, W0) per Sirlin 2011 Eq. (11)
    (Use same numerical stability measures as for g-function)
    
    R_nu_i = 1 + (alpha/(2*pi)) * h(W_hat, W0)
```

### 8.3 Conversion Algorithm (Electron → Neutrino Spectrum)

```
Algorithm: ElectronToNeutrinoConversion(W_nu, W0, Z, f_e)
Input:  W_nu — neutrino energies
        W0   — endpoint
        Z    — nuclear charge
        f_e  — function returning electron spectrum values
Output: f_nu — neutrino spectrum values

For each W_nu in array:
    W_hat = W0 - W_nu
    
    # Evaluate electron spectrum at W_hat
    f_e_val = f_e(W_hat)
    
    # Compute difference correction
    Compute delta_diff = (alpha/(2*pi)) * [h(W_hat, W0) - g(W_hat, W0)]
    per Sirlin 2011 Eq. (20)
    
    f_nu_i = f_e_val * (1 + delta_diff)
```

---

## 9. Physical Interpretation and Validation

### 9.1 Dominant Terms

The leading \(O(\alpha)\) correction is the **universal Sirlin function** \(g(W, W_0)\). Its dominant contribution comes from the logarithm:

\[3\ln\left(\frac{m_p}{m_e}\right) \approx 15.4\]

multiplied by \(\alpha/2\pi \approx 1.16 \times 10^{-3}\), giving a typical ~1.8% correction. The energy-dependent terms modulate this by ~±1% across the spectrum.

### 9.2 Sign Conventions

- The correction factor \(R(W, W_0)\) **multiplies** the uncorrected spectrum.
- \(R > 1\) means the correction increases the transition probability at that energy.
- For the electron spectrum: the correction typically **decreases** the probability at high \(W\) (\(W/W_0 \gtrsim 0.44\)) and **increases** it at low \(W\) (Sirlin 2011, Sec. III).
- For the neutrino spectrum: the correction is **positive and monotonically increasing** with \(W_\nu\) (Sirlin 2011, Sec. II.G).

### 9.3 Consistency Checks

1. **Lifetime check:** Integrate the corrected and uncorrected spectra:

   \[\int_1^{W_0} \frac{d\Gamma_0}{dW} R(W, W_0) dW\]

   The relative correction to the lifetime should be finite (the logarithmic divergence at the endpoint integrable over the spectrum).

2. **Electron-neutrino consistency:** The lifetime calculated from the corrected electron spectrum must equal that from the corrected neutrino spectrum.

3. **Zero-charge limit:** For \(Z = 0\), \(R(W, W_0) \to 1 + \frac{\alpha}{2\pi}g(W, W_0)\).

4. **Low-energy limit:** As \(W \to 1\) (non-relativistic), the correction should approach the values tabulated in Sirlin 1967.

---

## 10. Summary of Required Special Functions

| Function | Definition | Numerical Source |
|----------|------------|------------------|
| \(\text{Li}_2(x)\) | \(\sum_{k=1}^\infty x^k/k^2\) | Standard library (e.g., `scipy.special.spence` in Python, `gsl_sf_dilog` in C) |
| \(\tanh^{-1}(x)\) | \(\frac{1}{2}\ln\frac{1+x}{1-x}\) | Built-in or \(\text{arctanh}(x)\) |
| \(\gamma_E\) | Euler-Mascheroni constant | 0.57721566490153286 |

---

## References

1. A. Sirlin, Phys. Rev. 164, 1767 (1967) — Universal correction function \(g(E, E_m, m)\) and proof of model independence.
2. A. Sirlin and R. Zucchini, Phys. Rev. Lett. 57, 1994 (1986) — Analytic \(O(Z\alpha^2)\) calculation resolving CVC discrepancy.
3. A. Sirlin, Phys. Rev. D 35, 3420 (1987) — Refined finite-nuclear-size effects with realistic charge distributions.
4. A. Sirlin, Phys. Rev. D 84, 014021 (2011) — Neutrino spectrum radiative correction and conversion formula.
5. T. Kinoshita and A. Sirlin, Phys. Rev. 113, 1652 (1959) — Mass singularity cancellation theorems.
