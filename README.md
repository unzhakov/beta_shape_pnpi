# Beta Spectrum Calculation Toolkit

## Overview

This project provides a modular and extensible framework for calculating high-precision **beta decay energy spectra**. It is designed to model the spectral shape of a given nuclear transition by combining well-established physical corrections into a coherent and configurable pipeline.

The primary goal is to enable accurate numerical evaluation of the beta spectrum:

* from near-zero kinetic energies up to the endpoint
* with controlled inclusion of higher-order corrections
* using a clean, composable architecture

The framework focuses purely on **spectrum formation**, independent of any specific experimental or phenomenological application.

---

## Core Concept

The beta spectrum is constructed as a product of independent correction factors:

[
N(W) \propto p W (W_0 - W)^2 \times \prod_i C_i(W)
]

Each correction is implemented as a standalone component:

* phase space
* Coulomb interaction (Fermi function)
* finite nuclear size
* screening
* exchange
* radiative corrections
* recoil effects

This modular approach allows:

* easy validation of individual corrections
* selective enabling/disabling of physics effects
* flexible extension with new components

---

## Project Structure

```
beta_corr_pnpi/
├── beta_spectrum/
│   ├── base.py
│   ├── spectrum.py
│   ├── constants.py
│   ├── utils.py
│   │
│   ├── components/
│   │   ├── phase_space.py
│   │   ├── fermi.py
│   │   ├── finite_size.py
│   │   ├── screening.py
│   │   ├── exchange.py
│   │   ├── radiative.py
│   │   └── recoil.py
│   │
│   ├── data/
│   │   └── exchange_coeff.csv
│   │
│   └── __init__.py
│
├── tests/
│
├── pyproject.toml
├── README.md
```

---

## Architecture

### `SpectrumComponent` (base class)

Defined in `base.py`.

All physics corrections inherit from this class and implement:

```python
def __call__(self, W: np.ndarray) -> np.ndarray
```

This makes every component:

* stateless (or minimally stateful)
* vectorized
* composable

---

### `Spectrum`

Defined in `spectrum.py`.

Responsible for:

* combining all components
* evaluating the total spectrum
* managing energy grids

Typical usage pattern:

```python
spectrum = Spectrum(
    components=[
        PhaseSpace(W0),
        FermiFunction(Z),
        ScreeningCorrection(Z),
        ExchangeCorrection(Z),
        RadiativeCorrection(W0),
    ]
)

values = spectrum(W)
```

---

## Components

Each correction is implemented as an independent module under:

```
beta_spectrum/components/
```

---

### Phase Space (`phase_space.py`)

Implements:

[
p W (W_0 - W)^2
]

Provides the baseline spectral shape.

---

### Fermi Function (`fermi.py`)

Accounts for Coulomb interaction between the emitted electron and nucleus.

---

### Finite Size (`finite_size.py`)

Corrects for non-pointlike nuclear charge distribution.

---

### Screening (`screening.py`)

Includes atomic electron screening effects.

* modifies low-energy behavior
* becomes negligible at high energies

---

### Exchange (`exchange.py`)

Accounts for electron exchange with atomic orbitals.

* dominant at very low energies
* requires empirical coefficients (loaded from CSV)

---

### Radiative (`radiative.py`)

Implements outer radiative corrections:

* based on analytical expressions
* includes optional endpoint resummation
* handles logarithmic divergences near endpoint

---

### Recoil (`recoil.py`)

Includes nuclear recoil corrections.

---

## Data Files

### `exchange_coeff.csv`

Contains tabulated coefficients for exchange correction.

Loaded dynamically based on nuclear charge ( Z ).

---

## Design Principles

### 1. Modularity

Each physical effect is isolated in its own component.

---

### 2. Numerical Stability

Special care is taken near:

* ( W \to 1 ) (threshold)
* ( W \to W_0 ) (endpoint)

Typical strategies include:

* small epsilon regularization
* safe evaluation of logarithms and powers
* masking invalid regions

---

### 3. Vectorization

All components operate on NumPy arrays:

* no Python loops in hot paths
* efficient evaluation over dense grids

---

### 4. Extensibility

New corrections can be added by subclassing:

```python
class NewCorrection(SpectrumComponent):
    def __call__(self, W):
        return ...
```

---

## Typical Workflow

1. Define endpoint energy ( W_0 )
2. Construct energy grid
3. Initialize components
4. Combine into spectrum
5. Evaluate and plot

---

## Future Extensions

Planned improvements include:

* support for multiple decay branches
* nuclear shape factor implementations
* improved atomic corrections
* validation against reference spectra
* export utilities and CLI tools

---

## Notes

* Units are typically expressed in natural units ( m_e = 1 )
* Energies can be converted from MeV externally
* The framework assumes familiarity with beta decay theory

---

## License / Usage

(To be defined)

---

## Author Notes

This project is under active development and serves both as:

* a computational tool
* a structured reference for beta spectrum calculations

Updates may refine both physics implementations and code structure.
