# Beta Spectrum Calculation Toolkit

A modular Python package for calculating high-precision **beta decay energy spectra**. It models the spectral shape of nuclear transitions by combining well-established physical corrections into a coherent, configurable pipeline.

**Features:**

- Numerical evaluation from near-zero kinetic energies up to the endpoint
- Selective inclusion of higher-order corrections via `SpectrumConfig` toggles
- Clean, composable component architecture
- Built-in visualization and CSV export via `BetaSpectrumAnalyzer`
- Comprehensive scientific reference in [`docs/`](docs/)

______________________________________________________________________

## Core Concept

The beta spectrum is constructed as a product of independent correction factors:

$$N(W) \\propto p W (W_0 - W)^2 \\times \\prod_i C_i(W)$$

Each $C_i(W)$ is an independently implemented and toggleable component.

______________________________________________________________________

## Project Structure

```
beta_shape_pnpi/
├── pyproject.toml
├── docs/                          # Obsidian knowledge base (physics theory & references)
│   ├── README.md                  # Vault overview + navigation index
│   ├── 00-beta-spectrum-overview.md    # Master equation, all corrections table
│   ├── 01-fermi-function.md           # Fermi function F₀(Z,W)
│   ├── 02-phase-space.md              # Phase space p·W·(W₀−W)²
│   ├── 03-finite-size.md              # L₀, U, DFS corrections
│   ├── 04-screening-correction.md     # Atomic screening S(Z,W)
│   ├── 05-exchange-correction.md      # Electron exchange X(Z,W)
│   ├── 06-radiative-corrections.md    # Outer radiative δᵣ(W,W₀)
│   ├── 07-recoil-effects.md           # Weak magnetism, recoil kinematics
│   ├── 08-detector-response.md        # Detector response models
│   ├── 10-nuclear-structure.md        # Shape factor C(Z,W), BB vs HS formalisms
│   ├── 12-atomic-overlap.md           # Bahcall correction r(Z,W)
│   ├── 13-chemical-effects.md         # Molecular environment effects
├── beta_spectrum/
│   ├── __init__.py          # Public API: constants, utilities, classes
│   ├── base.py              # Abstract SpectrumComponent base class
│   ├── constants.py         # Physical constants (natural units, m_e = 1)
│   ├── utils.py             # Helpers: T<->W conversion, momentum, nuclear radius
│   ├── spectrum.py          # BetaSpectrum + BetaSpectrumAnalyzer classes
│   ├── nuclear_data.py      # paceENSDF integration + JSON input support
│   ├── cli.py               # Command-line interface (bs_pnpi)
│   ├── fitter.py            # χ² curve fitting for C(W) extraction
│   ├── cw_extractor.py      # C(W) shape factor extraction + gV/gA analysis
│   └── components/
│       │   ├── phase_space.py           ✓ Phase space shape (p·W·(W₀−W)²)
│       │   ├── fermi.py                 ✓ Coulomb correction (loggamma for stability)
│       │   ├── finite_size.py           ✓ L0 expansion + charge distribution U term
│       │   ├── screening.py             ✓ Atomic electron screening (ratio method)
│       │   ├── exchange.py              ✓ Hayen-2018 empirical fit coefficients
│       │   ├── radiative.py             ✓ Outer radiative corrections, soft-photon resummation
│       │   └── detector_response.py     ✓ Analytical detector response models
├── tests/
│   ├── test_radiative.py        # Radiative correction tests (24 tests)
│   ├── test_exchange.py         # Exchange correction tests
│   ├── test_fermi.py            # Fermi function tests
│   ├── test_finite_size.py      # Finite size correction tests
│   ├── test_phase_space.py      # Phase space tests
│   ├── test_screening.py        # Screening correction tests
│   ├── test_spectrum.py         # Integration tests
│   ├── test_detector_response.py# Detector response tests
│   ├── test_fitter.py           # Curve fitter tests
│   ├── test_cw_extractor.py     # C(W) extraction tests
│   └── test_nuclear_data.py     # paceENSDF + JSON input tests (23 tests)
└── data/
    ├── exchange_coeff.csv         # Tabulated coefficients for X(Z,W), Z=2..120
    └── custom_input_example.json  # Sample JSON input file
```

______________________________________________________________________

## Architecture

### `SpectrumComponent` (base class)

Defined in `beta_spectrum/base.py`. All physics corrections inherit from this abstract base and implement:

```python
def __call__(self, W: np.ndarray) -> np.ndarray
```

This makes every component stateless (or minimally stateful), vectorized, and composable.

### `BetaSpectrum`

Defined in `beta_spectrum/spectrum.py`. Orchestrates the calculation:

- Combines components multiplicatively via `__call__(W)`
- Factory method `from_config(config)` for declarative setup
- Generates energy grids with `get_energy_grid(config)`
- Returns individual component values via `calculate_components(W)`

### `BetaSpectrumAnalyzer`

Also in `spectrum.py`. Provides debugging and visualization tools:

- `total_spectrum(normalize=True/False)` — computed spectrum array
- `plot_analysis(save_path=None)` — 4-panel matplotlib figure (spectrum, components, cumulative effect, deviations from unity)
- `export_to_csv(filename)` — export via pandas DataFrame
- `get_data()` — return all numerical data for custom analysis

______________________________________________________________________

## Quick Start

```python
from beta_spectrum import SpectrumConfig, BetaSpectrum, BetaSpectrumAnalyzer

config = SpectrumConfig(
    Z_parent=90,       # Thorium-232 (example)
    Z_daughter=91,
    A_number=232,
    endpoint_MeV=4.8,  # Q-value
    transition_type="A",  # Allowed transition (default)
)

spectrum = BetaSpectrum.from_config(config)
W, E_MeV = spectrum.get_energy_grid(config)
values = spectrum(W)

analyzer = BetaSpectrumAnalyzer(spectrum, config)
analyzer.plot_analysis()          # Show/ save 4-panel figure
analyzer.export_to_csv("output.csv")   # Export data to CSV
```

______________________________________________________________________

## Components

| Component             | Module           | Status | Description                                                                                                    |
| --------------------- | ---------------- | -----: | -------------------------------------------------------------------------------------------------------------- |
| Phase Space           | `phase_space.py` |      ✓ | Baseline $p W (W_0 - W)^2$ with transition-type-dependent forbidden factors and optional neutrino mass support |
| Fermi Function        | `fermi.py`       |      ✓ | Relativistic Coulomb correction via `scipy.special.loggamma`                                                   |
| Finite Size L0        | `finite_size.py` |      ✓ | Low-Z expansion from Hayen et al. with prefactor                                                               |
| Charge Distribution U | `finite_size.py` |      ✓ | Second-order $(1/5)(\alpha Z W R)^2$ correction                                                               |
| Screening             | `screening.py`   |      ✓ | Ratio method with logistic switching, energy shift $V_0$                                                       |
| Exchange              | `exchange.py`    |      ✓ | Hayen-2018 Table X empirical fit (CSV data for Z=2..120)                                                       |
| Radiative             | `radiative.py`   |      ✓ | Outer radiative corrections (Hayen Eq. 47–53), soft-photon resummation at endpoint                             |
| Detector Response     | `detector_response.py` | ✓ | Analytical models: Gaussian, Gaussian+tail, Tikhonov, tabulated                                               |

## Additional Modules

| Module                | Description                                                                                                    |
| --------------------- | -------------------------------------------------------------------------------------------------------------- |
| `nuclear_data.py`     | paceENSDF integration + JSON input for automated `SpectrumConfig` generation                                   |
| `cli.py`              | Command-line interface (`bs_pnpi`) for quick spectrum calculations                                             |
| `fitter.py`           | χ² curve fitting for C(W) shape factor extraction from experimental data                                       |
| `cw_extractor.py`     | C(W) extraction pipeline with gV/gA coupling constant analysis                                                 |
