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

## 3. Test-Driven Development

### TDD Cycle

1. **Red** — write a failing test that expresses the desired behavior
2. **Green** — implement minimal code to make the test pass
3. **Refactor** — clean up code while keeping all tests green

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

## 6. Release Checklist

- [ ] All tests pass (`pytest`)
- [ ] Code is formatted (`black --check`)
- [ ] No lint errors (`ruff check`)
- [ ] Type checks pass (`mypy`)
- [ ] Update version in `pyproject.toml`
- [ ] Update `Development Status` section in `README.md`
- [ ] Commit changes
- [ ] Create git tag: `git tag v<version>`
- [ ] Push to remote
