# conftest.py — shared fixtures for pytest

import numpy as np
import pytest


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
    # Endpoint at T = 2.5 MeV → W0 ≈ 5.9, so we go up to W=5.8
    return np.linspace(1.003, 5.8, 100)


@pytest.fixture()
def large_W():
    """Large array for performance/shape checks."""
    return np.linspace(1.002, 4.9, 1000)
