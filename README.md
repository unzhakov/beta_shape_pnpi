# Beta Spectrum Calculation Toolkit

A modular Python package for calculating high-precision **beta decay energy spectra**. It models the spectral shape of nuclear transitions by combining well-established physical corrections into a coherent, configurable pipeline.

**Features:**
* Numerical evaluation from near-zero kinetic energies up to the endpoint
* Selective inclusion of higher-order corrections via `SpectrumConfig` toggles
* Clean, composable component architecture
* Built-in visualization and CSV export via `BetaSpectrumAnalyzer`
* Comprehensive scientific reference in [`docs/`](docs/)

---

## Core Concept

The beta spectrum is constructed as a product of independent correction factors:

$$N(W) \propto p W (W_0 - W)^2 \times \prod_i C_i(W)$$

Each $C_i(W)$ is an independently implemented and toggleable component.

---

## Project Structure

```
beta_shape_pnpi/
├── pyproject.toml
├── docs/                          # Obsidian knowledge base (physics theory & references)
│   ├── README.md                  # Vault overview + navigation index
│   ├── Hayen2017_summary.md       # Full paper summary backbone
│   ├── 00-beta-spectrum-overview.md    # Master equation, all corrections table
│   ├── 01-fermi-function.md           # Fermi function F₀(Z,W)
│   ├── 02-phase-space.md              # Phase space p·W·(W₀−W)²
│   ├── 03-finite-size.md              # L₀, U, DFS corrections
│   ├── 04-screening-correction.md     # Atomic screening S(Z,W)
│   ├── 05-exchange-correction.md      # Electron exchange X(Z,W)
│   ├── 06-radiative-corrections.md    # Outer radiative δᵣ(W,W₀)
│   ├── 07-recoil-effects.md           # Weak magnetism, recoil kinematics
│   ├── 10-nuclear-structure.md        # Shape factor C(Z,W), BB vs HS formalisms
│   ├── 12-atomic-overlap.md           # Bahcall correction r(Z,W)
│   ├── 13-chemical-effects.md         # Molecular environment effects
│   └── refs-hayen2017.md              # Key equations reference table, citations
├── beta_spectrum/
│   ├── __init__.py          # Public API: constants, utilities, classes
│   ├── base.py              # Abstract SpectrumComponent base class
│   ├── constants.py         # Physical constants (natural units, m_e = 1)
│   ├── utils.py             # Helpers: T<->W conversion, momentum, nuclear radius
│   ├── spectrum.py          # BetaSpectrum + BetaSpectrumAnalyzer classes
│   └── components/
│       │   ├── phase_space.py       ✓ Phase space shape (p·W·(W₀−W)²)
│       │   ├── fermi.py             ✓ Coulomb correction (loggamma for stability)
│       │   ├── finite_size.py       ✓ L0 expansion + charge distribution U term
│       │   ├── screening.py         ✓ Atomic electron screening (ratio method)
│       │   ├── exchange.py          ✓ Hayen-2018 empirical fit coefficients
│       │   ├── radiative.py         ✓ Outer radiative corrections, soft-photon resummation
│       │   └── recoil.py            ✗ Not yet implemented
└── data/
    └── exchange_coeff.csv     # Tabulated coefficients for X(Z,W), Z=2..120
```

---

## Architecture

### `SpectrumComponent` (base class)

Defined in `beta_spectrum/base.py`. All physics corrections inherit from this abstract base and implement:

```python
def __call__(self, W: np.ndarray) -> np.ndarray
```

This makes every component stateless (or minimally stateful), vectorized, and composable.

### `BetaSpectrum`

Defined in `beta_spectrum/spectrum.py`. Orchestrates the calculation:
* Combines components multiplicatively via `__call__(W)`
* Factory method `from_config(config)` for declarative setup
* Generates energy grids with `get_energy_grid(config)`
* Returns individual component values via `calculate_components(W)`

### `BetaSpectrumAnalyzer`

Also in `spectrum.py`. Provides debugging and visualization tools:
* `total_spectrum(normalize=True/False)` — computed spectrum array
* `plot_analysis(save_path=None)` — 4-panel matplotlib figure (spectrum, components, cumulative effect, deviations from unity)
* `export_to_csv(filename)` — export via pandas DataFrame
* `get_data()` — return all numerical data for custom analysis

---

## Quick Start

```python
from beta_spectrum import SpectrumConfig, BetaSpectrum, BetaSpectrumAnalyzer

config = SpectrumConfig(
    Z_parent=90,       # Thorium-232 (example)
    Z_daughter=91,
    A_number=232,
    endpoint_MeV=4.8,  # Q-value
)

spectrum = BetaSpectrum.from_config(config)
W, E_MeV = spectrum.get_energy_grid(config)
values = spectrum(W)

analyzer = BetaSpectrumAnalyzer(spectrum, config)
analyzer.plot_analysis()          # Show/ save 4-panel figure
analyzer.export_to_csv("output.csv")   # Export data to CSV
```

---

## Components

| Component | Module | Status | Description |
|---|---|---:|---|
| Phase Space | `phase_space.py` | ✓ | Baseline $p W (W_0 - W)^2$, optional neutrino mass support |
| Fermi Function | `fermi.py` | ✓ | Relativistic Coulomb correction via `scipy.special.loggamma` |
| Finite Size L0 | `finite_size.py` | ✓ | Low-Z expansion from Hayen et al. with prefactor |
| Charge Distribution U | `finite_size.py` | ✓ | Second-order $(1/5)(\alpha Z W R)^2$ correction |
| Screening | `screening.py` | ✓ | Ratio method with logistic switching, energy shift $V_0$ |
| Exchange | `exchange.py` | ✓ | Hayen-2018 Table X empirical fit (CSV data for Z=2..120) |
| Radiative | `radiative.py` | ✓ | Outer radiative corrections (Hayen Eq. 47–53), soft-photon resummation at endpoint |
| Recoil | `recoil.py` | ✗ | Stub only, not yet implemented |

---

## Configuration

Use `SpectrumConfig` to selectively enable/disable each correction:

```python
config = SpectrumConfig(
    Z_parent=20, Z_daughter=21, A_number=40, endpoint_MeV=5.0,
    use_fermi=True,
    use_screening=False,   # disable atomic screening
    use_exchange=True,
    use_radiative=True,
)
```

---

## Dependencies

**Runtime:** numpy, scipy (for `loggamma`, `spence`), matplotlib (plotting), pandas (CSV export)  
**Dev:** pytest, black, ruff, mypy (strict mode configured in `pyproject.toml`)

---

## Design Principles

### 1. Modularity
Each physical effect is isolated in its own component under `beta_spectrum/components/`.

### 2. Numerical Stability
Care near $W \to 1$ (threshold) and $W \to W_0$ (endpoint):
* `np.clip` / epsilon floors for regularization
* Safe logarithms via `scipy.special.loggamma` and dilogarithm
* Masking invalid regions

### 3. Vectorization
All components operate on NumPy arrays — no Python loops in hot paths.

### 4. Extensibility
New corrections are added by subclassing:

```python
class NewCorrection(SpectrumComponent):
    def __call__(self, W):
        return ...
```

---

## Development Status

**Version:** 0.1.2  
**Status:** Active development — core physics engine is functional; recoil correction and tests remain incomplete.

### Remaining Work

* Implement nuclear recoil correction in `recoil.py` (see [[07-recoil-effects]])
* Add unit / integration tests (pytest configured but no test files)
* Declare pandas as a runtime dependency in `pyproject.toml`
* Satisfy mypy strict type-checking requirements
* CLI entry point and multiple decay branch support

---

## Knowledge Base

The [`docs/`](docs/) directory contains the complete scientific reference for beta spectrum shape corrections, organized as an Obsidian vault. It is structured around Hayen et al., *Rev. Mod. Phys.* **90**, 015008 (2017) and maps directly to calculator components:

| # | Note | Calculator Module | Status |
|---|------|-------------------|--------|
| 1 | [[01-fermi-function]] | `components/fermi.py` | ✓ implemented |
| 2 | [[03-finite-size]] | `components/finite_size.py` | ✓ implemented |
| 3 | [[04-screening-correction]] | `components/screening.py` | ✓ implemented |
| 4 | [[05-exchange-correction]] | `components/exchange.py` | ✓ implemented |
| 5 | [[06-radiative-corrections]] | `components/radiative.py` | ✓ implemented |
| 6 | [[07-recoil-effects]] | `components/recoil.py` | ✗ not yet implemented |

Notes prefixed with `10-`, `12-`, `13-` cover nuclear structure, atomic overlap, and chemical effects — important for future extensions but outside the current calculator scope.

---

## Notes

* Physical constants are defined in natural units ($m_e = 1$); energies can be converted from MeV via utilities in `utils.py`
* The framework assumes familiarity with beta decay theory
