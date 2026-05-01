# conftest.py — shared fixtures for notebooks
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI

import matplotlib.pyplot as plt
import numpy as np
import os
import pytest

# Auto-save figures directory
PLOTS_DIR = os.path.join(os.path.dirname(__file__), "_plots")
os.makedirs(PLOTS_DIR, exist_ok=True)
_plot_counter = 0


def _save_figure_on_close(fig):
    """Save figure to PNG before closing."""
    global _plot_counter
    plot_path = os.path.join(PLOTS_DIR, f"notebook_{_plot_counter:03d}.png")
    try:
        fig.savefig(plot_path, dpi=100, bbox_inches="tight")
        _plot_counter += 1
    except Exception:
        pass


# Patch plt.close to auto-save before closing
_original_close = plt.close


def _patched_close(fig=None):
    if fig is not None:
        _save_figure_on_close(fig)
    return _original_close(fig)


plt.close = _patched_close
plt.rcParams.update({"figure.autolayout": True})


@pytest.fixture(autouse=True)
def close_figures():
    """Close all matplotlib figures after each notebook cell execution."""
    yield
    plt.close("all")


@pytest.fixture()
def W_low():
    """Low-energy total energy: T ≈ 10 keV → W = 1 + 0.01/0.511 ≈ 1.02."""
    return np.array([1.02])


@pytest.fixture()
def W_mid():
    """Mid-range energies (bulk of spectrum)."""
    return np.linspace(1.1, 2.0, 5)


@pytest.fixture()
def W_high():
    """High-energy total energy near endpoint: T ≈ 3 MeV → W ≈ 7."""
    return np.array([6.8])


@pytest.fixture()
def W_full():
    """Full spectrum grid from ~1 keV to just below endpoint."""
    return np.linspace(1.003, 5.8, 100)


@pytest.fixture()
def large_W():
    """Large array for performance/shape checks."""
    return np.linspace(1.002, 4.9, 1000)
