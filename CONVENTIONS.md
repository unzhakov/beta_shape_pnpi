# Project Conventions

This file defines the physical, mathematical, and architectural conventions that govern code correctness and consistency across the project.

## 1. Natural Units System

All internal representations, calculations, and definitions use natural units:

```
m_e = ℏ = c = 1
```

where `m_e` is the electron rest mass. This is the fundamental unit for energy, mass, momentum, and time throughout the codebase.

## 2. Energy Representation

### Core Principle

Internal calculations use natural units (m_e). External interfaces (API signatures, plots, human-readable output) use eV-based units (MeV or keV).

### API Signatures (Input)

Use human-readable MeV for parameters that come from databases or literature:

- `endpoint_MeV` — Q-value / endpoint energy in MeV
- `detector_sigma_a_keV` — resolution in keV
- `detector_tau_keV` — tail decay constant in keV

### Internal Representation

All internal computations use **total energy W in m_e units**:

```
W = T_to_W(T_MeV)   # T (kinetic) in MeV → W (total) in m_e
W = T + 1.0          # W = T/m_e + 1  (since m_e = 1 in m_e units)
```

- `W = 1.0` corresponds to zero kinetic energy (rest mass)
- `W > 1.0` corresponds to kinetic energy
- All internal functions accept and operate on `W`
- Endpoint W0 is always `T_to_W(endpoint_MeV)`

### Output and Display

Convert back to eV units for anything intended for human perception:

- Plots: x-axis in keV or MeV (kinetic energy)
- ASCII output, CSV export: kinetic energy in keV or MeV
- Log messages, error messages: use keV/MeV with context

```python
kinetic_keV = (W - 1.0) * ME_KEV    # W in m_e → kinetic energy in keV
kinetic_MeV = (W - 1.0) * ME_ME     # W in m_e → kinetic energy in MeV
```

### Conversion Helpers

All conversions and constants must go through the central utilities:

```python
from beta_spectrum.utils import T_to_W, W_to_T
from beta_spectrum.constants import ME_MEV, ME_KEV, ALPHA, MP_MEV
```

- Never hardcode physical constants — import from `constants.py`
- Never use manual arithmetic for unit conversion — use `T_to_W()`, `W_to_T()`
- Never mix MeV and m_e units in the same calculation

## 3. Energy Unit Convention

### MeV Throughout

All energy parameters in the public API use **MeV** as the standard unit:

- `endpoint_MeV` — always in MeV
- `e_step_MeV` — always in MeV
- `detector_sigma_a_keV` — in keV (detector-specific, exception noted)
- `detector_tau_keV` — in keV (detector-specific, exception noted)

Detector parameters use keV because detector resolutions are typically specified in keV in the literature. All other energies use MeV.

### CLI Unit Convention

All CLI energy parameters use MeV unless explicitly marked otherwise:

- `--e-step` — MeV (default: 0.001)
- `--sigma` — keV (detector resolution, default: 1.0)
- `--tau` — keV (tail decay constant, default: 5.0)

### Consistency Rule

When reading or writing values, always check the unit suffix in the parameter name:
- `_MeV` suffix → value is in MeV
- `_keV` suffix → value is in keV
- No suffix on energy parameter → value is in m_e units (internal)

## 4. SpectrumComponent Design

All spectral components (Fermi function, screening, radiative correction, etc.) follow a common interface:

- **Stateless or minimally stateful** — no mutable state that depends on input
- **Vectorized** — accept and return `np.ndarray`, support array broadcasting
- **Consistent energy variable** — always use `W` (total energy in m_e units) internally
- **Callable interface** — implement `__call__(W)` for uniform usage
- **Return correction factor** — output is a multiplicative factor, typically ~0.1 to 10

```python
class MyComponent(SpectrumComponent):
    def __init__(self, param1: float, param2: float):
        self.param1 = param1
        self.param2 = param2

    def __call__(self, W: np.ndarray) -> np.ndarray:
        # All internal computation uses W in m_e units
        return ...
```

## 5. Logging Convention

### Verbosity Levels

- **No verbosity (default)**: Only errors and warnings
- **`-v` (INFO)**: Shows source selection, component initialization, workflow steps, and output locations
- **`-vv` (DEBUG)**: Shows all INFO messages plus detailed internals of each spectrum component (parameters, input/output ranges, evaluation counts)

### Component Logging

Each `SpectrumComponent` logs at two levels:

- **INFO** at initialization: Component name and key parameters
- **DEBUG** at evaluation: Number of energy points, input/output value ranges

This allows debugging without cluttering normal output.

## 6. Nuclear Data Notation

### Element and Nuclide Notation

- Element symbol: `Tc`, `Ru`, `U`, etc. (from atomic number)
- Nuclide notation: `Tc99`, `Ru99`, `U238` (element symbol + mass number)
- Decay notation: `Tc99 -> Ru99` (parent -> daughter)

### CSV Header Convention

CSV metadata headers use element notation for readability:

```
# nuclide: Tc99 -> Ru99 (Z=43->44)
```

This provides both the human-readable notation and the machine-parseable Z values.

## 7. Physics-Specific Testing Guidelines

When writing tests for physics components, verify the following:

- **Physical constraints** — values at thresholds, endpoints, and boundaries must match known physics behavior
- **Monotonicity / positivity** — spectrum components should not produce unphysical negative values; monotonicity should match physical expectations
- **Shape correctness** — output shape must match input shape (vectorized operation)
- **Type safety** — components accept and return `np.ndarray`
- **Numerical stability** — no NaN or inf in the physical energy range, except where physically expected (e.g., below kinematic threshold, near soft-photon endpoint divergence)
