Here is the updated summary with the inline citations integrated as requested.

# Summary: The measurement of the $^{99}$Tc $\beta$-decay spectrum and its implications for the effective value of weak axial coupling

## 1. Introduction and Motivation

The study of $\beta$-spectral shapes is a critical method for determining the effective value of the weak axial coupling constant, $g_A$. In forbidden non-unique $\beta^-$ transitions [5], the spectral shape is directly sensitive to the ratio $g_A/g_V$.

* **Standard Hypothesis:** The Conserved Vector Current (CVC) hypothesis sets the vector coupling constant $g_V = 1.0$ [8].
* **Quenching:** The Partially Conserved Axial-Vector Current (PCAC) hypothesis suggests that in finite nuclei, $g_A$ is "quenched" (reduced) to a value of approximately $1.0$, compared to the bare-nucleon value of $1.27$ resulting from the experimental analysis of the decay of an isolated neutron [8]. The quenching of the effective value of $g_A$ has been studied extensively in [2, 8, 9], and in particular its effects on the neutrinoless double beta decay have been addressed in [1, 3].
* **Contradiction:** While previous studies on the fourth-forbidden non-unique $\beta^-$-decays of $^{113}$Cd [10–13] and $^{115}$In [14–16] supported quenching ($g_A/g_V \approx 1.0$), a recent measurement by Paulsen et al. [17] on the second-forbidden non-unique transition in $^{99}$Tc suggested an *enhanced* ratio of $g_A/g_V = 1.526(92)$.

This paper aims to resolve this contradiction by performing a high-precision measurement of the $^{99}$Tc $\beta$-decay spectrum using Metallic Magnetic Calorimeters (MMCs) and comparing it with state-of-the-art nuclear theory.

## 2. Theoretical Framework

The half-life ($t_{1/2}$) and shape of the $\beta$-spectrum are governed by the integrated shape function $\tilde{C}$.

**Equation 1:**
$$ t_{1/2} = \frac{\kappa}{\tilde{C}} $$
*(where $\kappa$ is a constant [18, 19] and $\tilde{C}$ incorporates phase-space factors and NMEs arising in the next-to-leading-order expansion, as discussed in detail in Refs. [6, 7].)*

The shape function $C(w_e)$ depends on the weak axial and vector couplings:

**Equation 2:**
$$ C(w_e) = g_V^2 C_V(w_e) + g_A^2 C_A(w_e) + g_V g_A C_{VA}(w_e) $$
*(where $w_e$ is the total energy of the emitted electron. The transition matrix element [6].)*

**Methodology (ESSM):**
The authors employ the **Enhanced Spectrum-Shape Method (ESSM)** [19]. This approach simultaneously fits the measured $\beta$-spectrum and the experimental partial half-life to extract $g_A^{eff}$. This allows for the adjustment of the small relativistic vector nuclear matrix element (sNME) [5] to reproduce the half-life, rather than fixing it to the "ideal" CVC value.

The sNME gathers contributions outside the valence major shell, making its calculation particularly hard for nuclear-theory frameworks used, e.g., in the present work and also in the work of Paulsen et al. [17]. Despite its smallness, the sNME can influence the $\beta$-spectral shapes and half-lives quite strongly [11, 12, 19, 21].

A realistic value of the sNME can be estimated from its CVC value [5], first used for spectral-shape calculations in Ref. [18] and later by Paulsen et al. [17]. To circumvent the use of the CVC value, the experimental works [11, 13, 15] and theoretical works [12, 16] used sNME as a fitting parameter, together with $g_A^{eff}$.

In this Letter, we present new theoretical calculations of the $\beta$ spectrum of $^{99}$Tc... In references [21, 23] the related two sets of theoretical $\beta$-spectra of $^{99}$Tc were introduced for two NSM Hamiltonians, namely $jj45pnb$ [24] and $glekpn$ [25].

## 3. Experimental Setup

The experiment was conducted at Lawrence Livermore National Laboratory (LLNL).

* **Detector:** A Metallic Magnetic Calorimeter (MMC) fabricated at KRISS [26, 27].
  * **Sensor:** Paramagnetic ions (420 ppm of Ag:$^{168}$Er) embedded in a metallic host.
  * **Absorber:** A 4$\pi$ solid-angle gold absorber containing the $^{99}$Tc source.
  * **Temperature:** Operated below 100 mK (specifically 15 mK during data taking).
* **Source Preparation:** $^{99}$Tc solution was dropped onto a gold foil, dried, folded, and rolled (kneading method [28, 29]) to ensure uniform distribution. The final sample mass was 1.75 mg.
* **Readout:** A Magnicon SQUID (Model X1) and an additional array SQUID (Model X16).
* **Data Acquisition:** Continuous mode, 200 kS/s sampling rate using a Stanford Research Systems SRS 560 amplifier [30].
* **Threshold:** A low-energy threshold of 10 keV was applied.

**Figure 1: Detector Setup**
*(Description of visual content)*
The image shows the experimental apparatus. Key components labeled include:

* **Magnicon SQUID X1:** The signal readout device.
* **KRISS MMC:** The calorimeter sensor.
* **$^{99}$Tc source:** The radioactive sample placed on the sensor.
* **Support washer & Kapton tape:** Used for thermal insulation and alignment.

## 4. Results and Analysis

### 4.1 Data Selection

* **Background Dataset:** 130,011 events recorded.
* **Signal Dataset:** After Pulse Shape Discrimination (PSD) cuts, **77,487 events** corresponding to $^{99}$Tc remained for analysis.

### 4.2 Energy Calibration

Calibration was performed using an external $^{133}$Ba source and internal X-rays from a $^{148}$Gd source (used for internal calibration checks).

**Table I: Measured X-ray and gamma-ray lines**

| Origin | Radiation type | Tabulated energy [keV] | Measured energy [keV] | FWHM resolution [keV] | $\Delta$E [keV] |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Gd fluorescence** | X-ray K$_{\alpha2}$ | 42.3097(13) | 42.30(10) | 0.9(4) | 0.01(10) |
| | X-ray K$_{\alpha1}$ | 42.9968(12) | 42.9427(31) | 0.68(10) | 0.05(33) |
| | X-ray K$_{\beta3}$ | 48.5505(21) | 48.29(5)* | 1.05(22) | 0.41(5) |
| | X-ray K$_{\beta1}$ | 48.6957(21) | | | |
| **Pt fluorescence** | X-ray K$_{\alpha2}$ | 65.1226(21) | 65.005(22) | 0.9(6) | 0.118(22) |
| | X-ray K$_{\alpha1}$ | 66.8311(20) | 66.78(5) | 0.65(15) | 0.06(5) |
| **$^{133}$Ba calibration** | $\gamma$ | 80.9979(11) | 80.881(24) | 0.77(8) | 0.117(25) |
| | $\gamma$ | 276.3989(12) | 276.12(5) | 0.57(13) | 0.28(5) |
| | $\gamma$ | 302.8508(5) | 302.764(21) | 0.60(5) | 0.087(21) |
| | $\gamma$ | 356.0129(7) | 355.926(9) | 0.653(21) | 0.087(9) |
| | $\gamma$ | 383.8485(12) | 383.69(6) | 0.82(14) | 0.16(6) |

*(Note: Entries labeled "Not distinguished" in the original text refer to peaks with low relative intensities that could not be individually resolved.)*

### 4.3 Spectral Comparison

The measured spectrum was compared against theoretical predictions from three nuclear models:

1. **$jj45pnb$** (NSM Hamiltonian)
2. **$glekpn$** (NSM Hamiltonian)
3. **MQPM** (Microscopic Quasiparticle-Phonon Model)

**Figure 2: $^{99}$Tc spectrum with external $^{133}$Ba calibration source**
*(Description)*
A plot showing "Counts / 200 keV" vs "Energy (keV)". It displays the raw data (black line) with distinct peaks corresponding to Ba-133 $\gamma$-lines, Au X-rays, Pt X-rays, and Gd X-rays [34].

**Figure 3: Spectrum of calibrated $^{99}$Tc compared with theoretical predictions**
*(Description)*
A plot showing "Counts / 1 keV" vs "Energy (keV)" (0-300 keV).

* **Black line:** Experimental Data.
* **Blue dashed line:** $jj45pnb$ model ($g_A=1.0$).
* **Green dotted line:** $glekpn$ model ($g_A=1.2$).
* **Red dash-dot line:** MQPM model ($g_A=1.0$).
* **Bottom panel:** Residuals ($\Delta/c$) showing the fit quality over the 0-280 keV range.

### 4.4 Quantitative Results (Table II)

The table below summarizes the fit results for the effective axial coupling $g_A^{eff}$ and the goodness of fit ($\chi^2_\nu$ and p-value).

**Table II: Results for the three considered nuclear models**

| Model | $g_A^{eff}$ | sNME [fm$^3$] | $\chi^2_\nu$ | p-value |
| :--- | :--- | :--- | :--- | :--- |
| **CloseCVC** | | | | |
| glekpn | 1.1 | 0.0674 | 1.1057 | 0.12 |
| jj45pnb | 1.0 | 0.0681 | 1.2436 | 0.004 |
| MQPM | 1.0 | 0.0651 | 1.6437 | <0.0001 |
| **FarCVC (Best fit)** | | | | |
| glekpn | 1.2 | -0.0698 | 1.0664 | 0.22 |
| jj45pnb | 1.0 | -0.0669 | 1.0966 | 0.14 |
| MQPM | 1.0 | -0.0683 | 1.0951 | 0.14 |

* **CloseCVC:** Assumes the sNME is fixed to the CVC value. The $jj45pnb$ and MQPM models show poor agreement (low p-values).
* **FarCVC (Best fit):** Allows sNME to vary. The **$glekpn$ model with $g_A^{eff} = 1.2$** yields the best agreement ($\chi^2_\nu = 1.0664$, p-value = 0.22). The other models with $g_A^{eff} = 1.0$ also show statistically compatible fits (p-values $\approx$ 0.14).

## 5. Conclusion

The high-precision measurement of the $^{99}$Tc $\beta$-decay spectrum supports a **quenched axial coupling** ($g_A^{eff} \approx 1.0 - 1.2$).

* This contradicts the previously reported enhanced value ($g_A/g_V = 1.526$) by Paulsen et al. [17].
* The discrepancy is attributed to the different treatment of the sNME (small relativistic vector nuclear matrix element) in the theoretical analysis.
* The results are consistent with quenching patterns observed in heavier nuclei, such as $^{113}$Cd [10–13] and $^{115}$In [14–16].

## References

1. J. Engel and J. Menendez, Rep. Prog. Phys. **60**, 046301 (2017).
2. H. Ejiri, J. Suhonen, and K. Zuber, Phys. Rep. **797**, 1 (2019).
3. M. Agostini, G. Benato, J. A. Detwiler, J. Menéndez, and F. Vissani, Rev. Mod. Phys. **95**, 025002 (2023).
4. J. Suhonen, *From Nucleons to Nucleus: Concepts of Microscopic Nuclear Theory* (Springer, Berlin, 2007).
5. H. Behrens and W. Bühring, *Electron Radial Wave Functions and Nuclear Beta-decay* (International Series of Monographs on Physics) (Clarendon Press, Oxford, 1982).
6. M. Haaranen, P. C. Srivastava, and J. Suhonen, Phys. Rev. C **93**, 034308 (2016).
7. M. Haaranen, J. Kotila, and J. Suhonen, Phys. Rev. C **95**, 024327 (2017).
8. J. Suhonen, Front. Phys. **5**, 55 (2017).
9. J. Suhonen and J. Kotela, Front. Phys. **7**, 29 (2019).
10. L. Bodenstein-Dresler *et al.* (The COBRA Collaboration), Phys. Lett. B **800**, 135092 (2020).
11. J. Kostensalo, J. Suhonen, J. Volkmer, S. Zatschler, and K. Zuber, Phys. Lett. B **822**, 136652 (2021).
12. J. Kostensalo, E. Lisi, A. Marrone, and J. Suhonen, Phys. Rev. C **107**, 055502 (2023).
13. I. Bandac *et al.*, Eur. Phys. J. C **(2024)** 84:1158.
14. A. F. Leder *et al.*, Phys. Rev. Lett. **129**, 232502 (2022).
15. L. Pagnanini *et al.* (The ACCESS Collaboration), Phys. Rev. Lett. **133**, 122501 (2024).
16. J. Kostensalo, E. Lisi, A. Marrone, and J. Suhonen, Phys. Rev. C **110**, 055503 (2024).
17. M. Paulsen *et al.*, Phys. Rev. C **110**, 055503 (2024).
18. A. Kumar, P. C. Srivastava, J. Kostensalo, and J. Suhonen, Phys. Rev. C **101**, 064304 (2020).
19. A. Kumar, P. C. Srivastava, and J. Suhonen, Eur. Phys. J. A **57**, 225 (2021).
20. O. Nitescu, S. Stoica, and F. Šimkovic, Phys. Rev. C **107** (2023) 025501.
21. M. Ramalho and J. Suhonen, Phys. Rev. C **109**, 034321 (2024).
22. M. Ramalho, J. Suhonen, A. Neacsu, and S. Stoica, Front. Phys. **12**:1455778 (2024).
23. M. Ramalho and J. Suhonen, Il Nuovo Cimento **47** C, 377 (2024).
24. A. F. Lisetskiy, B. A. Brown, M. Horoi, and H. Grawe, Phys. Rev. C **70**, 044314 (2004).
25. H. Mach *et al.*, Phys. Rev. C **41**, 226 (1990).
26. J. W. Song, S. G. Kim, H. S. Kim, H. J. Kim, and M. K. Lee, J. Low Temp. Phys. **216**, 436 (2024).
27. J. W. Song, Y. C. Cho, H. J. Kim, and M. K. Lee, J. Low Temp. Phys. **218**, 110 (2025).
28. M. P. Croce, A. S. Hoover, M. W. Rabin, E. M. Bond, L. E. Wolfsberg, D. R. Schmidt, and J. N. Ullom, J. Low Temp. Phys. **184**, 938 (2016).
29. A. S. Hoover *et al.*, Anal. Chem. **87**, 3996 (2015).
30. Stanford Research Systems, SR560 Low-Noise Preamplifier, Rev. 2.8 (Stanford Research Systems, Sunnyvale, CA, 2006), <https://www.thinksrs.com/products/sr560.html>.
31. V. T. Jordanov, G. F. Knoll, A. C. Huber, and J. A. Pantazis, Nucl. Instrum. Methods Phys. Res., Sect. A **353**, 261 (1994).
32. A. R. L. Kavner, D. Lee, S. T. P. Boyd, S. Friedrich, I. Jovanovic, and G. B. Kim, J. Low Temp. Phys. **209**, 1070 (2022).
33. The MathWorks, Inc., MATLAB, version R2016b (The MathWorks, Inc., Natick, MA, 2016).
34. R. D. Deslattes, E. G. Kessler Jr., P. Indelicato, L. De Billy, E. Lindroth, and J. Anton, Rev. Mod. Phys. **75**, 35 (2003).
35. Y. Khazov, A. Rodionov, and F. G. Kondev, Nucl. Data Sheets **112**, 855 (2011).
36. A. C. Thompson *et al.*, X-Ray Data Booklet (Center for X-Ray Optics and Advanced Light Source, Lawrence Berkeley National Laboratory, Berkeley, CA, 2009).
