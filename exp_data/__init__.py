"""
Experimental data handling for beta spectrometry.

Provides tools for loading, calibrating, and pre-processing
experimental spectra from detector readouts.
"""

from .spectrum import ExpSpectrum
from .calibration import EnergyCalibrator
from .fitters import GaussianFitter, PeakFitter
from .corrections import (
    DeadTimeCorrection,
    PileUpCorrection,
    BackgroundSubtractor,
)

__all__ = [
    "ExpSpectrum",
    "EnergyCalibrator",
    "GaussianFitter",
    "PeakFitter",
    "DeadTimeCorrection",
    "PileUpCorrection",
    "BackgroundSubtractor",
]
