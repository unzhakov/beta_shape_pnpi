# Agent Instructions — beta-spectrum

High-precision beta decay spectrum calculator. Extracts shape factor C(W) from ⁹⁹Tc beta spectrum to parametrize g_V coupling.

## Stack & Commands

```
pytest                           # all tests
nbmake notebooks/                # execute notebooks as tests
black --check .                  # formatting
ruff check .                     # linting
mypy .                           # type checking (strict)
bs_pnpi --nuclide Tc99 -vv --output ./output/verification.csv  # debug verification
```

Install: `source ~/.pyenv/sci/bin/activate && pip install -e ".[dev]"`

## Git Workflow (GitHub Flow)

- **Never commit directly to `main`.** Always branch: `git checkout -b dev/<short-desc>`
- **Atomic commits** with **Conventional Commits** format:
  - `feat:`, `fix:`, `docs:`, `style:`, `refactor:`, `perf:`, `test:`, `chore:`
  - Imperative mood, present tense, under 72 chars
  - WIP commits fine on dev branches — squash before final merge
- **Commit each logical change** — don't wait for green. Use `git add -p`.
- **Rebase over merge** for linear history: `git rebase main` before merge.
- Delete branch locally + remotely after merge.

## TDD

Test-first by default. For new classes/functions: write failing test first.

- **RED** — failing test defining desired behavior
- **GREEN** — minimal code to pass
- **REFACTOR** — clean up, confirm tests still pass

State which phase you're in. Name phases explicitly.

### Python Testing Conventions

- Tests mirror package structure: `beta_spectrum/components/foo.py` → `tests/physics/test_foo.py`
- Shared fixtures in `tests/conftest.py`
- Classes organized by feature area (e.g., `TestComponentBasicProperties`)
- Descriptive names: `test_<function>_<expected_behavior>_<context>`
- Include implementation-level assertions (type checks, pre/postconditions)
- Use `hypothesis` for property-based tests where applicable
- Physics tests verify: physical constraints at thresholds, positivity, monotonicity, numerical stability

## Quality Gates

All must pass before merge:

```bash
pytest
nbmake notebooks/
black --check .
ruff check .
mypy .
```

Run individually to identify failures.

## Debug Verification (Final Step)

After quality gates pass, run:

```bash
bs_pnpi --nuclide Tc99 -vv --output ./output/verification.csv --plot ./output/verification.png --log-file ./output/debug.log 2>&1 | tee ./output/debug_stdout.log
```

Verify parameter consistency across all components. Save artifacts to `./output/`.
