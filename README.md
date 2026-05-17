# Beta Spectrum Calculation Toolkit

A modular Python package for calculating high-precision **beta decay energy spectra**. It models the spectral shape of nuclear transitions by combining well-established physical corrections into a coherent, configurable pipeline.

**Features:**

- Numerical evaluation from near-zero kinetic energies up to the endpoint
- Selective inclusion of higher-order corrections via `SpectrumConfig` toggles
- Clean, composable component architecture
- Built-in visualization and CSV export via `BetaSpectrumAnalyzer`
  - Normal mode: single spectrum plot with nuclear data header
  - Debug mode (`-vv`): 4-panel view with all correction factors
  - All plots include commit hash and UTC timestamp for traceability
- Comprehensive scientific reference in [`docs/`](docs/)
- Debug verification workflow: `bs_pnpi -vv` validates all parameters end-to-end

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
├── beta_spectrum/           ← pure theory core
│   ├── __init__.py          # Public API: constants, utilities, classes
│   ├── base.py              # Abstract SpectrumComponent base class
│   ├── constants.py         # Physical constants (natural units, m_e = 1)
│   ├── utils.py             # Helpers: T<->W conversion, momentum, nuclear radius
│   ├── spectrum.py          # BetaSpectrum + SpectrumConfig + BranchConfig
│   ├── nuclear_data.py      # paceENSDF integration + JSON input support
│   ├── cli.py               # Command-line interface (bs_pnpi)
│   ├── logging_utils.py     # Structured logging helpers
│   └── components/
│       │   ├── phase_space.py           ✓ Phase space shape (p·W·(W₀−W)²)
│       │   ├── fermi.py                 ✓ Coulomb correction (loggamma for stability)
│       │   ├── finite_size.py           ✓ L0 expansion + charge distribution U term
│       │   ├── screening.py             ✓ Atomic electron screening (ratio method)
│       │   ├── exchange.py              ✓ Hayen-2018 empirical fit coefficients
│       │   ├── radiative.py             ✓ Outer radiative corrections, soft-photon resummation
│       │   └── detector_response.py     ✓ Analytical detector response models
├── beta_spectrum/visualize/ ← plotting / visualization
│   └── __init__.py          # BetaSpectrumAnalyzer
├── exp_data/                ← experimental data handling
│   ├── __init__.py          # ExpSpectrum, EnergyCalibrator, GaussianFitter, etc.
│   ├── spectrum.py          # ExpSpectrum dataclass (energies, counts, errors, metadata)
│   ├── calibration.py       # EnergyCalibrator (linear/quadratic channel→energy)
│   ├── fitters.py           # GaussianFitter + PeakFitter for cal sources
│   └── corrections.py       # DeadTimeCorrection, PileUpCorrection, BackgroundSubtractor
├── fitter/                  ← analysis framework
│   ├── __init__.py          # SpectrumModel, SpectrumFitter, CWExtractor, etc.
│   ├── model.py             # SpectrumModel (BetaSpectrum + DetectorResponse wrapper)
│   ├── fit_engine.py        # SpectrumFitter + CurveFitter + FitConfig + FitResult
│   ├── extractor.py         # CWExtractor, GVAExtractor, CWExtractionResult, GVAExtractionResult
│   └── result.py            # AnalysisResult + FitSummary
├── tests/
│   ├── exp_data/            # exp_data module tests
│   ├── fitter/              # fitter module tests
│   ├── physics/             # Physics component tests
│   ├── quality/             # Code quality tests (logging, API, nuclear data)
│   ├── integration/         # Full pipeline integration tests
│   └── common/              # Shared hypothesis + property tests
├── data/
│   ├── exchange_coeff.csv         # Tabulated coefficients for X(Z,W), Z=2..120
│   └── custom_input_example.json  # Sample JSON input file
├── output/                      # Test run outputs (gitignored except .gitkeep)
├── AGENTS.md                    # Agent instructions: stack, git workflow, TDD, quality gates
├── CONVENTIONS.md               # Physics conventions: natural units, energy representation, component design
├── TODO.md                      # Development roadmap
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

In `beta_spectrum.visualize`. Provides debugging and visualization tools:

- `total_spectrum(normalize=True/False)` — computed spectrum array
- `plot_analysis(save_path=None, show_components=True)` — plot with two modes:
  - `show_components=True` (debug mode): 4-panel figure (spectrum, components, cumulative effect, deviations)
  - `show_components=False` (normal mode): single spectrum plot with nuclear data header
  - All plots include commit hash and UTC timestamp for traceability
- `export_to_csv(filename)` — export via pandas DataFrame with YAML-style metadata header
- `get_data()` — return all numerical data for custom analysis

______________________________________________________________________

## Quick Start

```python
from beta_spectrum import SpectrumConfig, BetaSpectrum
from beta_spectrum.visualize import BetaSpectrumAnalyzer

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
analyzer.plot_analysis("analysis.png")          # Normal mode (single spectrum plot)
analyzer.plot_analysis("debug.png", show_components=True)  # Debug mode (4-panel)
analyzer.export_to_csv("spectrum.csv")   # Export data to CSV
```

### Fitting and C(W) extraction

```python
from beta_spectrum import SpectrumConfig, BetaSpectrum
from fitter import SpectrumModel, SpectrumFitter, CWExtractor
from exp_data.spectrum import ExpSpectrum

# Build theoretical model
spectrum = BetaSpectrum.from_config(config)
model = SpectrumModel(spectrum, detector_response=detector)

# Fit to experimental data
exp = ExpSpectrum(energies=energies_keV, counts=counts, errors=errors)
fitter = SpectrumFitter(model, exp.energies, exp.counts, exp.errors)
result = fitter.fit(x0=[1.0, 0.0])

# Extract C(W) shape factor
extractor = CWExtractor(config)
cw_result = extractor.extract(measured, errors, energies_keV, model_values, endpoint_keV=Q)
```

______________________________________________________________________

## Implementation Status

For a complete list of implemented and planned features, see [TODO.md](TODO.md). Key highlights:

- **Physics corrections**: phase space, Fermi function, finite size, screening, exchange, radiative (with delta_cut resummation)
- **Detector response**: Gaussian, Gaussian+tail, Tikhonov models with tabulated support
- **Multi-branch decay**: intensity-weighted sum of branch spectra, per-branch calculators with branch-specific endpoint, transition type (from ENSDF), and `--intensity-cutoff` CLI filter
- **Analysis tools**: χ² fitter (SpectrumFitter + CurveFitter), C(W) extraction pipeline, Kurie plot analysis, g_V/g_A extraction
- **Experimental data**: ExpSpectrum container, energy calibration, peak fitting, dead-time/pile-up correction
- **CLI**: `bs_pnpi --nuclide Tc99` — full paceENSDF integration (auto-lookup of Z, A, Q-value, transition type, half-life, branches, forbiddenness from ENSDF)
- **paceENSDF integration**: `get_decay_info_from_paceENSDF()`, `create_config_from_source("paceENSDF", nuclide="Tc99")`, `_parse_nuclide_symbol()` for Z=1..98, FORBIDDENNESS_MAP, `--nuclide` CLI flag
