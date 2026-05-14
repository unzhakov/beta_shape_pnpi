# First-forbidden transitions in the reactor anomaly

**Authors:** L. Hayen, J. Kostensalo, N. Severijns, and J. Suhonen
**Affiliations:**

1. Instituut voor Kern- en Stralingsfysica, KU Leuven, Belgium
2. Department of Physics, University of Jyväskylä, Finland
**Date:** August 23, 2019
**Source:** arXiv:1908.08302v1 [nucl-th]

---

## Abstract

This paper describes microscopic calculations of dominant forbidden transitions in reactor antineutrino spectra above 4 MeV using the nuclear shell model. By incorporating Coulomb corrections in the most complete way, the authors calculate the shape factor with high fidelity, showing strong deviations from allowed approximations. Despite small differences in ab initio electron cumulative spectra, large differences (several percent) are found in antineutrino spectra. The authors propose a parametrization of forbidden spectra and derive spectral corrections and uncertainties using Monte Carlo techniques. They conclude that a correct treatment of forbidden transitions is indispensable for solving both the reactor normalization anomaly and the spectral shoulder.

## I. Introduction

The field of short baseline reactor neutrinos has seen significant activity, driven by long-standing issues like the LSND anomaly [1, 2], GALLIX & SAGE collaborations [3], and the Reactor Antineutrino Anomaly (RAA) [4, 5]. The RAA involves a normalization deficit and a "5 MeV bump" in the spectrum [8–12]. Theoretical predictions often rely on the Huber and Mueller et al. treatments [14, 15], which introduced strong approximations for forbidden transitions. This work aims to provide a more complete analysis using the nuclear shell model for transitions above 4 MeV.

## II. $\beta$ Decay Formalism

The treatment of forbidden $\beta$ decays is complex, involving kinematic, nuclear, and Coulomb terms. The transition matrix element is given by:

$$
\mathcal{M}_{fi} = \int \mathrm{d}^3r \, \bar{\phi}_e(\vec{r}, \vec{p}_e) \gamma^\mu (1+\gamma^5) v(\vec{p}_\nu) \times \int \frac{\mathrm{d}^3s}{(2\pi)^3} e^{i\vec{s}\cdot\vec{r}} \frac{1}{2} [\langle f | \hat{p}(f) + \vec{p}_e - \vec{s} | V_\mu + A_\mu | i \rangle \langle \vec{p}_i \rangle + \langle f | \hat{p}(f) | V_\mu + A_\mu | i \rangle \langle \vec{p}_e + \vec{s} \rangle]
$$

The $\beta$ spectrum shape is traditionally written as:

$$
\frac{\mathrm{d}N}{\mathrm{d}W} = \frac{G_V^2 V_{ud}^2}{2\pi^3} p W (W-W_0)^2 \times F(Z,W) C(Z,W) K(Z,W)
$$

where $C(Z,W)$ is the shape factor. The shape factor depends on:

1. Spin-change and kinematic factors.
2. Finite size corrections proportional to $R^n$.
3. Coulomb corrections proportional to $(\alpha Z)^n$.

### A. Nuclear structure

Forbidden transitions correspond to $\beta$ decays where Fermi and Gamow-Teller matrix elements are zero. The authors use the Behrens-Bühring formalism [23] to label nuclear structure form factors using quantum numbers $K, L, s$. For first-forbidden $\beta$ transitions ($\Delta J = 0, 1, 2$ and parity change $\pi_i \pi_f = -1$), up to 6 form factors contribute.

### B. Coulomb corrections

Coulomb corrections are split into static Coulomb renormalization and convolution distortion. The authors note that neglecting electron mass and Coulomb interaction leads to a symmetric shape factor, an argument often used to neglect forbidden transitions in the RAA context, which they argue is invalid.

### C. Breakdown of usual approximations

The $\xi$ approximation, often used to simplify shape factors, is shown to be questionable for relevant fission fragments where $\alpha Z / 2R \sim W_0$.

## III. Data Selection & Handling

The study relies on the ENDF/B-VIII.0 decay data library [38] and ENSDF data. To address the "Pandemonium effect" [37], where high-energy $\gamma$ decays are missed, the authors use TAGS (Total Absorption Gamma Spectroscopy) data [34–36] to correct branching ratios.

## IV. Shape Factor Calculation

The authors calculated shape factors for 36 dominant forbidden transitions using the NUSHELLX@MSU shell model code [48].

**Table I: Dominant forbidden transitions above 4 MeV**
*(Selected entries from Table I)*

| Nuclide | $Q_\beta$ (MeV) | $E_{ex}$ (MeV) | BR (%) | $J_i^\pi \to J_f^\pi$ | FY (%) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $^{89}$Br | 8.3 | 0 | 16 | $3/2^- \to 3/2^+$ | 1.1 |
| $^{90}$Rb | 6.6 | 0 | 33 | $0^- \to 0^+$ | 4.5 |
| $^{91}$Kr | 6.8 | 0.11 | 18 | $5/2^+ \to 5/2^-$ | 3.5 |
| $^{92}$Rb | 8.1 | 0 | 95.2 | $0^- \to 0^+$ | 4.8 |
| ... | ... | ... | ... | ... | ... |

The calculated shape factors (Fig. 3) show significant deviations from unity, especially for pseudovector transitions.

## V. Spectral Changes

Comparisons with existing literature (Fang and Brown [17]) show good agreement for $^{140}$Xe but discrepancies for $^{136}$Te, highlighting the sensitivity to nuclear structure details. The inclusion of forbidden shape factors leads to a downward trend in the electron spectrum below 4 MeV and an increase at higher energies (Fig. 5).

## VI. Improved Forbidden Transition Treatment

To generalize results to the entire database, the authors propose a parametrization of the shape factor:

$$
C = 1 + aW + b/W + cW^2
$$

They use Gaussian kernel density estimation to determine the distributions of fit parameters $a, b, c$ (Fig. 8, 9) and employ Monte Carlo techniques to estimate uncertainties.

## VII. Updated Summation Calculations

The authors combine their results with the ILL data set [47]. They compare three approaches:

1. Treating transitions as allowed with weak magnetism correction.
2. Treating them as allowed with $C=1$.
3. Using the calculated forbidden shape factors.

**Figure 14:** Relative change to cumulative antineutrino spectra. A bump appears between 4 and 7 MeV with a magnitude of up to 4.5%.

**Table III:** Difference in integral and IBD flux compared to Huber-Mueller results.

| | $^{235}$U | $^{238}$U | $^{239}$Pu | $^{241}$Pu |
| :--- | :--- | :--- | :--- | :--- |
| $\phi$ | 0.2(2) | 0.4(5) | 0.2(2) | 0.3(2) |
| $R_{IBD}$ | 0.8(5) | 2.3(10) | 0.7(5) | 0.7(6) |

## VIII. Reactor Spectrum Changes

The study investigates the spectral shoulder (Fig. 17) and rate anomaly (Fig. 18). The authors find that while forbidden transitions alone do not fully solve the anomaly, combining their results with a proposed increase in the average slope of allowed shape factors can solve both the rate anomaly and the spectral shoulder simultaneously.

## IX. Conclusion

The authors conclude that:

1. Forbidden transitions above 4 MeV are dominant and show strong deviations from allowed approximations.
2. A correct treatment of these transitions is essential for accurate reactor antineutrino predictions.
3. The parametrization of forbidden shape factors provides a robust method for uncertainty estimation.
4. The combination of forbidden transition effects and slope changes in allowed transitions offers a potential solution to the reactor anomaly and spectral bump.

---

## References

1. C. Athanassopoulos et al., *Physical Review Letters* **81**, 1774 (1998).
2. J. M. Conrad, W. C. Louis, and M. H. Shaevitz, *Annual Review of Nuclear and Particle Science* **63**, 45 (2013).
3. F. Kaether et al., *Physics Letters B* **785**, 570 (2018).
4. A. A. Aguilar-Arevalo et al., *Physical Review D* **94**, 112008 (2016).
5. G. Mention et al., *Physical Review D* **83**, 073006 (2011).
6. K. N. Abazajian et al., *arXiv:1204.5309* (2012).
7. J. R. Jordan et al., *Physical Review Letters* **122**, 081801 (2019).
8. A. C. Hayes et al., *Physical Review D* **92**, 033015 (2015).
9. P. Huber, *Physical Review Letters* **118**, 042502 (2017).
10. F. P. An et al., *Physical Review Letters* **116**, 061801 (2016).
11. Y. Abe et al., *Physical Review Letters* **116**, 061801 (2016).
12. S. H. Seo et al., *arXiv:1610.04326* (2016).
13. A. C. Hayes and P. Vogel, *Annual Review of Nuclear and Particle Science* **66**, 219 (2016).
14. P. Huber, *Physical Review C* **84**, 024611 (2011); *Physical Review C* **85**, 029901 (2012).
15. T. A. Mueller et al., *Physical Review C* **83**, 055515 (2011).
16. A. C. Hayes et al., *Physical Review Letters* **112**, 202501 (2014).
17. D.-L. Fang and B. A. Brown, *Physical Review C* **91**, 025503 (2015).
18. L. Hayen et al., *Physical Review C* **99**, 031301(R) (2019).
19. B. Stech and L. Schülke, *Zeitschrift für Physik* **179**, 314 (1964).
20. B. R. Holstein, *Physical Review C* **19**, 1544 (1979).
21. L. Hayen et al., *Reviews of Modern Physics* **90**, 015008 (2018).
22. H. Behrens and W. Bühring, *Nuclear Physics A* **162**, 111 (1971).
23. H. Behrens and W. Bühring, *Electron radial wave functions and nuclear beta-decay* (Clarendon Press, Oxford, 1982).
24. E. K. Warburton and I. S. Towner, *Physics Reports* **242**, 103 (1994).
25. P. Baumann et al., *Physical Review C* **58**, 1970 (1998).
26. H. Behrens and J. Jänecke, *Landolt-Börnstein Tables, Gruppe I, Band 4* (Springer, 1969).
27. M. Morita and R. S. Morita, *Physical Review* **109**, 2048 (1958).
28. H. F. Schopper, *Weak Interactions and Nuclear Beta Decay* (North-Holland Publishing Company, 1966).
29. T. Kotani, *Physical Review* **114**, 795 (1959).
30. W. Bühring, *Nuclear Physics A* **430**, 1 (1984).
31. B. Holstein, *Reviews of Modern Physics* **46**, 789 (1974).
32. X. B. Wang and A. C. Hayes, *Physical Review C* **95**, 064313 (2017).
33. M. Fallot et al., *arXiv:1208.3877* (2012).
34. A. Algara et al., *Physical Review Letters* **105**, 202501 (2010).
35. A.-A. Zakari-Issoufou et al., *Physical Review Letters* **115**, 102503 (2015).
36. S. Rice et al., *arXiv:1702.05512* (2017).
37. J. C. Hardy et al., *Physics Letters B* **711**, 307 (1977).
38. M. Navratil et al., *Nuclear Data Sheets* **148**, 1 (2018).
39. S. C. Van de Marck, *Nuclear Data Sheets* **113**, 2935 (2012).
40. M. Kellett, O. Bersillon, and R. Mills, *The JEFF-3.1-3.1.1 RadioActive Decay Data and Fission Yields Sub-Libraries* (JEFF Report 20, 2009).
41. M. Fallot and A. Sonzogni, *Private Communication*.
42. M. Herman and A. Trkov, *Brookhaven National Laboratory, Tech. Rep.* (2010).
43. M. Chadwick et al., *Nuclear Data Sheets* **112**, 2887 (2011).
44. T. Kawano, P. Möller, and W. B. Wilson, *Physical Review C* **78**, 054601 (2008).
45. A. A. Sonzogni, T. D. Johnson, and E. A. McCutchan, *Physical Review C* **91**, 011301 (2015).
46. L. Hayen and N. Severijns, *Computer Physics Communications* **240**, 152 (2019).
47. N. Haag et al., *arXiv:1405.3554* (2014).
48. B. A. Brown and W. D. M. Rae, *Nuclear Data Sheets* **120**, 115 (2014).
49. H. Mach et al., *Physical Review C* **41**, 226 (1990).
50. R. Machleidt, *Physical Review C* **63**, 024001 (2001).
51. S. Lalkovski et al., *arXiv:1212.4961* (2012).
52. B. A. Brown, *Unpublished* (2012).
53. J. Suhonen, *Frontiers in Physics* **5**, 55 (2017).
54. J. Kostensalo and J. Suhonen, *Physics Letters B* **781**, 480 (2018).
55. J. Kostensalo, M. Haaranen, and J. Suhonen, *Physical Review C* **95**, 044313 (2017).
56. M. Haaranen, J. Kotila, and J. Suhonen, *Physical Review C* **95**, 024327 (2017).
57. L. Bodenstedt-Dresler et al., *arXiv:1806.02254* (2018).
58. J. P. Davidson, *Collective Models of the Nucleus* (Academic Press, New York and London, 1968).
59. J. Suhonen and J. Kostensalo, *Frontiers in Physics* **7**, 29 (2019).
60. A. A. Sonzogni, E. A. McCutchan, and A. C. Hayes, *Physical Review Letters* **119**, 112501 (2017).
61. D.-L. Fang, *Private Communication*.
62. T. Yoshida, T. Tachibana, S. Okumura, and S. Chiba, *Physical Review C* **98**, 41303 (2018).
63. D. W. Scott, *Multivariate Density Estimation: Theory, Practice, and Visualization* (John Wiley & Sons, 1992).
64. D. Foreman-Mackey, *The Journal of Open Source Software* **1**, 24 (2016).
65. D. A. Dwyer and T. J. Langford, *Physical Review Letters* **114**, 012502 (2015).
66. M. Estienne et al., *arXiv:1904.08358* (2019).
67. P. Vogel and J. F. Beacom, *Physical Review D* **60**, 10 (1999).
68. A. Strumia and F. Vissani, *Physics Letters B* **564**, 42 (2003).
69. A. Kurylov, M. J. Ramsey-Musolf, and P. Vogel, *Physical Review C* **67**, 035502 (2003).
70. F. P. An et al., *Physical Review Letters* **108**, 171803 (2012).
71. J. K. Ahn et al., *Physical Review Letters* **108**, 191802 (2012).
