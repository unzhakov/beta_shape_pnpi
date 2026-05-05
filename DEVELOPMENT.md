# Development Workflow

## 1. Setup

Do not attempt system wide pip installations. Use virtual environment 'sci'.
Default venv location: `~/.pyenv/`.

```bash
source ~/.pyenv/sci/bin/activate
pip install -e ".[dev]"
```

## 2. Git Workflow

### Before Starting Any Work

```bash
git status          # MUST show clean working tree
git fetch origin
git rebase origin/main   # update to latest main
```

If `git status` is not clean, commit or stash existing changes first. Never begin new work on a dirty working tree.

### Branch Strategy

```bash
git checkout main
git checkout -b dev/<feature-name>   # always create a new branch
```

Branch naming: `dev/<short-description>` (e.g., `dev/recoil-correction`)


### Merging Rules

- All work happens on `dev/*` branches — never commit directly to `main`
- Tests MUST pass before any merge attempt
- Merging to `main` requires a pull request with passing CI checks
- If tests fail, fix them on the dev branch — do NOT merge broken code
- Before merging always update TODO.md and README.md

## 3. Test-Driven Development

### TDD Cycle

1. **Red** — write a failing test that expresses the desired behavior
2. **Green** — implement minimal code to make the test pass
3. **Refactor** — clean up code while keeping all tests green

Commit work-in-progress changes to development branch.

### Test Conventions

- Tests live in `tests/` with the same structure as `beta_spectrum/`
- New component `beta_spectrum/components/foo.py` → test file `tests/test_foo.py`
- Use shared fixtures from `tests/conftest.py` (`W_low`, `W_mid`, `W_high`, `W_full`, `large_W`)
- Organize tests in classes by feature area (e.g., `TestComponentBasicProperties`, `TestComponentEdgeCases`)
- Each test class should have a docstring explaining what is being tested and why
- Use descriptive test names: `test_threshold_vanishes`, not `test_1`

### What to Test

See [`CONVENTIONS.md`](CONVENTIONS.md) §4 for physics-specific testing guidelines (physical constraints, positivity, numerical stability).

For general test structure, use the conventions below.

## 4. Quality Gates

All of the following MUST pass before merging to `main`:

```bash
pytest                           # all tests pass
nbmake notebooks/                # execute and validate notebooks
black --check .                  # code formatting
ruff check .                     # linting
mypy .                           # type checking (strict mode)
```

These are configured in `pyproject.toml`. Run them individually to identify failures.

### Notebook Quality Control

Notebooks in `notebooks/` are executed as pytest tests using `nbmake`. This ensures:

- All code cells execute without errors
- Full error traces are available for debugging
- Plots are extracted for LLM-based visual analysis

**Workflow:**

1. Write notebooks with explicit figure saving:
   ```python
   fig.savefig('_plots/notebook_name_cell{N}.png', dpi=100, bbox_inches='tight')
   from IPython.display import display, Image
   display(Image(filename='_plots/notebook_name_cell{N}.png'))
   ```

2. Execute and extract plots:
   ```bash
   python scripts/run_notebooks.py notebooks/
   ```

3. Send extracted plots to LLM for visual review:
   ```
   @notebooks/_plots/notebook_name_cell4.png
   ```

4. Run full quality gates:
   ```bash
   pytest notebooks/ --nbmake  # or just `pytest` (includes notebooks/)
   ```

______________________________________________________________________

### 4.5. Debug Verification (Final Step of Each Dev Iteration)

After all quality gates pass, run a **debug-level verification** as the final step of every development iteration. This catches configuration mismatches, missing parameter propagation, and silent logic errors before merging.

**Purpose:** Verify that all parameters (endpoint energy, energy grid, Z, A, physical constants) are consistent throughout the pipeline, and that every spectrum component emits sufficient debug logging for explicit comparison.

**Procedure:**

1. **Run `bs_pnpi` with `-vv` (DEBUG logging) and output to `./output/`:**
   ```bash
   bs_pnpi --nuclide Tc99 -vv --output ./output/verification.csv
   ```

2. **Capture full debug output to a log file for analysis:**
   ```bash
   bs_pnpi --nuclide Tc99 -vv --output ./output/verification.csv --log-file ./output/debug.log 2>&1 | tee ./output/debug_stdout.log
   ```

3. **Verify consistency of key parameters in the debug output:**

   | Parameter | What to check | Where to find it |
   |---|---|---|
   | `endpoint_MeV` | Matches expected Q-value | CLI dry-run, config creation, energy grid |
   | `Z_parent` / `Z_daughter` | Consistent across all components | Component init logs |
   | `A_number` | Consistent across all components | Component init logs |
   | `e_step_MeV` | Matches grid size × (endpoint - e_step) | Energy grid creation |
   | `W0` (endpoint in natural units) | `W0 = 1.0 + endpoint_MeV / m_e` | PhaseSpace init |
   | Physical constants | `m_e = 0.51099895` MeV, `α = 1/137.036` | Constants module, Fermi function |
   | Energy grid bounds | `e_step` to `endpoint - e_step` | `get_energy_grid()` output |

4. **Verify each component's debug logging:**

   Every enabled component must log its initialization parameters and a representative evaluation. Check the debug output for entries like:

   ```
   DEBUG ... Initializing PhaseSpace: W0=1.57664, transition_type=A, m_e=1.0
   DEBUG ... Evaluating PhaseSpace at 294 energy points
   DEBUG ... PhaseSpace range: [min=..., max=...]
   ```

   Required log entries per component:

   | Component | Must log |
   |---|---|
   | `PhaseSpace` | `W0`, `transition_type`, `m_e`, `m_nu`, energy range |
   | `FermiFunction` | `Z`, `A`, `α`, Coulomb approximation used |
   | `FiniteSizeL0` | `Z`, `A`, nuclear radius `R`, L0 approximation |
   | `ChargeDistributionU` | `Z`, `A`, charge radius, U correction range |
   | `ScreeningCorrection` | `Z_parent`, screening model, S(Z,W) range |
   | `ExchangeCorrection` | `Z_parent`, Hayen coefficients used, X(Z,W) range |
   | `RadiativeCorrection` | `W0`, `δ_r` formula, resummation flag, `delta_cut` value |

5. **Report any missing logging:**

   If a component's debug output does not allow explicit verification of its parameters or computed range, **add logging to that component** before the iteration is considered complete. The goal is that a developer (or an LLM) reading the debug log can reconstruct every parameter and verify correctness without reading source code.

6. **Save debug artifacts to `./output/`:**

   After verification, ensure all artifacts are in `./output/`:
   - `debug_stdout.log` — full stdout+stderr from `bs_pnpi -vv`
   - `debug.log` — log file output (if `--log-file` was used)
   - `verification.csv` — the spectrum CSV with metadata header
   - `verification.png` — analysis plot (if `--plot` was used)

   These files serve as a **verifiable record** of the last working state and can be used for regression comparison.

**Example verification command:**

```bash
# Full debug run with all outputs
bs_pnpi --nuclide Tc99 -vv \
    --output ./output/verification.csv \
    --plot ./output/verification.png \
    --log-file ./output/debug.log 2>&1 | tee ./output/debug_stdout.log

# Quick parameter extraction from debug log
grep -E '(Initializing|Evaluating|W0=|Z=|endpoint|e_step|PhaseSpace|Fermi|FiniteSize|Screening|Exchange|Radiative)' ./output/debug_stdout.log
```

**Success criteria:**

- [ ] All parameters (Z, A, endpoint, e_step, W0) appear consistently in debug output
- [ ] Each enabled component logs its initialization parameters
- [ ] Each enabled component logs its evaluation energy range
- [ ] Physical constants match expected values
- [ ] No silent parameter coercion or unit conversion errors detected
- [ ] All artifacts saved to `./output/`

**If any check fails:**

1. Add missing debug logging to the component
2. Fix any parameter mismatch
3. Re-run the verification until all checks pass
4. Only then consider the iteration complete

## 5. Agent Instructions

When working on this project, agents must:

1. **Read existing tests first** — understand patterns before writing new ones
2. **Follow existing conventions** — match test style, fixture usage, and class organization
3. **Never commit changes** without explicit user request
4. **Never merge branches** without explicit user request
5. **Verify `git status` is clean** before starting any file modifications
6. **Work on a dev branch** — create `dev/<description>` and switch to it
7. **Run quality gates** after making changes and report results
8. **Update docs** when changing public API or adding features
9. **Run debug verification** (`bs_pnpi -vv`) as the final step of each iteration, saving all artifacts to `./output/` and verifying parameter consistency across all components

## 6. Release Checklist

- [ ] All tests pass (`pytest`)
- [ ] Code is formatted (`black --check`)
- [ ] No lint errors (`ruff check`)
- [ ] Type checks pass (`mypy`)
- [ ] Debug verification passes (`bs_pnpi -vv`) — all parameters consistent, all components logged
- [ ] Debug artifacts saved to `./output/` (debug_stdout.log, verification.csv, etc.)
- [ ] Update version in `pyproject.toml`
- [ ] Update `Development Status` section in `README.md`
- [ ] Update `TODO.md` marking corresponding items as complete.
- [ ] Commit changes
- [ ] Create git tag: `git tag v<version>`
- [ ] Push to remote
