______________________________________________________________________

title: "08 – Detector Response Convolution"
date: 2026-05-01
tags:

- beta-spectrum/detector-response
  status: active
  aliases: [Convolution, Response Matrix, Smearing]
  cssclasses: []
  related_components: \[[00-beta-spectrum-overview]\], \[[06-radiative-corrections]\]

______________________________________________________________________

# Detector Response Convolution — $M(E_j) = \\sum_i \\Phi(W_i) \\cdot R(E_j; W_i) \\cdot \\Delta W_i$

## Overview

Detector response convolution bridges the gap between theoretical beta spectra and experimentally measured data. A theoretical spectrum $\\Phi(W)$ assumes perfect energy resolution, while real detectors smear the true energy distribution due to finite resolution, charge collection inefficiencies, and electronic noise.

The convolution operation:

$$M_j = \\sum_i \\Phi(W_i) \\cdot R(E_j; W_i) \\cdot \\Delta W_i$$

where:

- $\\Phi(W_i)$ — theoretical differential spectrum at true energy $W_i$
- $R(E_j; W_i)$ — detector response probability: probability that an electron of true energy $W_i$ deposits energy $E_j$ in channel $j$
- $\\Delta W_i$ — energy bin width at $W_i$
- $M_j$ — predicted measured counts in detector channel $j$

## Two Modes of Operation

### Tabulated Response

Receive a pre-computed response matrix from Monte Carlo simulations (GEANT4, EGSNRC) or dedicated calibration measurements. The matrix has shape $(N\_{\\text{calib}}, N\_{\\text{channels}})$, where each row $i$ is the detector response to mono-energetic electrons at calibration energy $E_i$.

```python
from beta_spectrum import DetectorResponse

# Load from MC simulation output
detector = DetectorResponse.from_mc_simulation(
    channel_energies=channel_grid,
    response_matrix=mc_response,
    calibration_energies=mc_energies,
)
```

### Analytical Model

Parameterized model with Gaussian core + exponential tail, suitable for quick estimates when no tabulated response is available.

```python
detector = DetectorResponse.from_gaussian_params(
    channel_energy_range=(1.0, 1.6),  # in m_e units
    n_channels=4096,
    sigma_a=0.0,       # energy-independent resolution (m_e)
    sigma_b=0.0,       # sqrt(E) contribution (m_e^{-1/2})
    tail_fraction=0.0, # fraction of events in exponential tail
    tau=0.01,          # decay constant of tail (m_e)
    model="gaussian",  # "gaussian" | "gaussian_tail" | "tikhonov"
    fano_factor=0.0,   # Fano factor (0.12 for Si, ~0 for MMC)
)
```

## Analytical Response Models

### Gaussian Core

Pure Gaussian peak with energy-dependent resolution:

$$G(E; E_0, \\sigma) = \\frac{1}{\\sigma\\sqrt{2\\pi}} \\exp\\left(-\\frac{(E - E_0)^2}{2\\sigma^2}\\right)$$

### Gaussian + Exponential Tail

Models low-energy tailing from charge collection inefficiencies in semiconductor detectors (standard for Si-PIPS detectors):

$$f(E) = (1 - f\_{\\text{tail}}) \\cdot G(E; E_0, \\sigma) + f\_{\\text{tail}} \\cdot T(E; E_0, \\sigma, \\tau)$$

where $T$ is the Tikhonov function:

$$T(E) = \\frac{1}{2\\tau} \\exp\\left(\\frac{\\tau(E - E_0)}{2} + \\frac{\\sigma^2\\tau^2}{8}\\right) \\cdot \\text{erfc}\\left(\\frac{E - E_0}{\\sqrt{2}\\sigma} + \\frac{\\sigma\\tau}{\\sqrt{2}}\\right)$$

### Energy-Dependent Resolution

$$\\sigma(E) = \\sqrt{\\sigma_a^2 + (\\sigma_b \\sqrt{E})^2 + F \\cdot E \\cdot w}$$

- $\\sigma_a$ — energy-independent contribution (electronic noise, etc.)
- $\\sigma_b$ — $\\sqrt{E}$ term (statistical fluctuations in charge creation)
- $F$ — Fano factor ($\\approx 0.12$ for Si, $\\approx 0$ for MMCs)
- $w$ — mean ionization energy ($\\approx 3.6$ eV for Si, converted to $m_e$ units)

## Calculator Implementation

> [!tip] Code mapping
> Module: `beta_spectrum/components/detector_response.py` → `DetectorResponse` class
>
> - **Analytical models**: Gaussian core, Gaussian + exponential tail (Tikhonov)
> - **Tabulated mode**: Linear interpolation from MC-calculated response matrix
> - **Convolution**: `DetectorResponse.convolve(W, spectrum, normalize=True)`
> - **Integration with BetaSpectrum**: `BetaSpectrum.convolve_with_detector(detector, W, config)`
> - **Integration with BetaSpectrumAnalyzer**: `BetaSpectrumAnalyzer.convolved_spectrum(detector, normalize)`

### API Usage

```python
from beta_spectrum import (
    BetaSpectrum, SpectrumConfig, DetectorResponse,
    BetaSpectrumAnalyzer, T_to_W
)

# 1. Create theoretical spectrum for 99Tc
config = SpectrumConfig(
    Z_parent=43, Z_daughter=44, A_number=99,
    endpoint_MeV=0.294  # 99Tc Q-value
)
spectrum = BetaSpectrum.from_config(config)
W, kinetic_MeV = spectrum.get_energy_grid(config)

# 2. Define detector response (sigma = 1.0 keV in m_e units)
sigma_1keV_me = 1.0 / 511.0  # convert keV to m_e units
W0 = T_to_W(0.294)
detector = DetectorResponse.from_gaussian_params(
    channel_energy_range=(1.0, W0 + 0.05),
    n_channels=512,
    sigma_a=sigma_1keV_me,
    sigma_b=0.0,
    tail_fraction=0.0,
    model="gaussian",
    fano_factor=0.0,
)

# 3. Convolve
convolved = spectrum.convolve_with_detector(detector, W=W, config=config)

# 4. Or use the Analyzer
analyzer = BetaSpectrumAnalyzer(spectrum, config)
convolved_normalized = analyzer.convolved_spectrum(detector, normalize=True)
```

## Convolution Algorithm

The convolution proceeds in three steps:

1. **Build response matrix**: For each energy bin $W_i$, compute the detector response $R(E_j; W_i)$ across all channels $j$. Result: matrix of shape $(N\_{\\text{bins}}, N\_{\\text{channels}})$.

1. **Weight by spectrum**: Multiply each response row by the theoretical spectrum value $\\Phi(W_i)$ and bin width $\\Delta W_i$.

1. **Sum over energies**: For each channel $j$, sum contributions from all energy bins:
   $$M_j = \\sum_i \\Phi(W_i) \\cdot R(E_j; W_i) \\cdot \\Delta W_i = (R^T \\cdot (\\Phi \\cdot \\Delta W))\_j$$

The bin width is computed as the half-bin average: $\\Delta W_i = \\frac{1}{2}(\\Delta W\_{i-1/2} + \\Delta W\_{i+1/2})$.

> [!note] Normalization
> Each response function $R(E; W_i)$ is normalized to unit area before convolution. This ensures that convolving a mono-energetic line produces a response with the same total counts. The final convolved spectrum preserves the integral of the theoretical spectrum.

## Energy Units

All detector response computations use **total energy in units of $m_e c^2$** (natural units, $m_e = 1$):

- $m_e c^2 = 510.998950$ keV
- Conversion: $W = 1 + T / (m_e c^2)$ where $T$ is kinetic energy
- A resolution of $\\sigma = 1.0$ keV corresponds to $\\sigma = 1.0 / 511.0 \\approx 0.00196$ in $m_e$ units

> [!warning] Unit consistency
> The detector channel energy range and the spectrum energy grid must be in the same units (total energy in $m_e$ units). Use `T_to_W()` from `beta_spectrum.utils` to convert kinetic energy in MeV to total energy in $m_e$ units.

## Physical Interpretation

The convolution operation smears the sharp theoretical spectrum into a realistic measured spectrum:

- **Peak broadening**: Finite energy resolution widens spectral features
- **Endpoint smearing**: The sharp cutoff at $W_0$ becomes a gradual fall-off
- **Tail effects**: Charge collection inefficiencies create low-energy tails (modeled by `gaussian_tail`)
- **Count conservation**: Total number of events is preserved (integral of convolved spectrum = integral of theoretical spectrum)

## References

- Hayen et al., RevModPhys **90**, 015008 (2017) — Comprehensive review of beta-decay theory
- Tikhonov, Soviet Physics Doklady **8**, 591 (1964) — Original formulation of the response function
- Sigmund, Phys. Rev. **130**, 1497 (1963) — Low-energy tailing in semiconductor detectors
