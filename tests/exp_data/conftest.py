"""Shared fixtures for exp_data tests."""

from __future__ import annotations

import numpy as np
import pytest

from exp_data.spectrum import ExpSpectrum
from exp_data.calibration import CalibrationResult


@pytest.fixture
def synthetic_spectrum():
    """Create a synthetic ExpSpectrum (4096 channels)."""
    n_channels = 4096
    energies = np.arange(n_channels, dtype=np.float64)
    counts = np.exp(-energies / 500) * 1000 + np.random.poisson(5, n_channels)
    return ExpSpectrum(
        energies=energies,
        counts=counts,
        source="Tc99",
        run_id="synthetic",
    )


@pytest.fixture
def calibration_result():
    """Create a synthetic CalibrationResult (linear, E = 0.4 + 0.1*ch)."""
    coeffs = np.array([0.1, 0.4])  # polyfit order: [slope, intercept]
    cov = np.eye(2) * 1e-6
    return CalibrationResult(
        coefficients=coeffs,
        order=1,
        chi2_per_dof=0.5,
        n_points=2,
        covariance=cov,
    )


@pytest.fixture
def full_spectrum():
    """Create a full ExpSpectrum with all fields populated."""
    return ExpSpectrum(
        energies=np.arange(4096, dtype=np.float64),
        counts=np.exp(-np.arange(4096) / 500) * 1000 + np.random.poisson(5, 4096),
        source="Tc99",
        run_id="run_042",
        date="2025-01-15",
        dead_time=0.05,
        live_time=3420.0,
        metadata={
            "calibration": {"type": "linear", "coefficients": [0.1, 0.4]},
            "livetime": {"fraction": 0.95, "method": "pulser"},
        },
    )
