# Project Conventions

Physical, mathematical, and architectural conventions governing code correctness.

## 1. Natural Units

All internal calculations use natural units:

```
m_e = ℏ = c = 1
```

`m_e` is the electron rest mass. Fundamental unit for energy, mass, momentum, and time.

## 2. Energy Representation

**Internal:** total energy `W` in m_e units (`W = 1.0` = rest mass, `W > 1.0` = kinetic).

**External API:** human-readable MeV/keV. Check the unit suffix in parameter names:
- `_MeV` suffix → value is in MeV
- `_keV` suffix → value is in keV (detector params only)
- No suffix on energy → value is in m_e units (internal)

### Conversion

Always use central utilities — never hardcode constants or manual arithmetic:

```python
from beta_spectrum.utils import T_to_W, W_to_T
from beta_spectrum.constants import ME_MEV, ME_KEV, ALPHA, MP_MEV
```

## 3. SpectrumComponent Design

All spectral components follow a common interface:

- **Stateless or minimally stateful** — no mutable state depending on input
- **Vectorized** — accept/return `np.ndarray`, support broadcasting
- **Consistent energy variable** — always `W` (total energy in m_e units) internally
- **Callable interface** — implement `__call__(W)` for uniform usage
- **Return correction factor** — multiplicative factor, typically ~0.1 to 10

## 4. Logging

- **No verbosity (default):** errors and warnings only
- **`-v` (INFO):** source selection, component init, workflow steps, output locations
- **`-vv` (DEBUG):** all INFO + detailed internals (parameters, ranges, evaluation counts)

Each `SpectrumComponent` logs:
- **INFO** at init: component name + key parameters
- **DEBUG** at evaluation: energy point count, input/output ranges

## 5. Nuclear Data Notation

- Element: `Tc`, `Ru`, `U`
- Nuclide: `Tc99`, `Ru99`, `U238`
- Decay: `Tc99 -> Ru99` (parent → daughter)
- CSV headers use element notation: `# nuclide: Tc99 -> Ru99 (Z=43->44)`

## 6. Physics Testing Guidelines

When testing physics components, verify:

- **Physical constraints** — values at thresholds, endpoints, boundaries match known physics
- **Positivity** — no unphysical negative values
- **Shape correctness** — output shape matches input (vectorized)
- **Type safety** — components accept/return `np.ndarray`
- **Numerical stability** — no NaN/inf in physical range (except where physically expected)
