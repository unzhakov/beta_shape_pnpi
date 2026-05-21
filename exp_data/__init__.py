"""
Experimental data handling for beta spectrometry.

Provides tools for loading, calibrating, and pre-processing
experimental spectra from detector readouts.
"""

from .spectrum import ExpSpectrum
from .calibration import EnergyCalibrator, CalibrationResult
from .fitters import GaussianFitter, PeakFitter
from .corrections import (
    DeadTimeCorrection,
    PileUpCorrection,
    BackgroundSubtractor,
)
from .livetime import LivetimeDetermination, LivetimeResult
from .pipeline import (
    SpectrumPipeline,
    PipelineResult,
    process_spectrum,
    process_run_sequence,
)
from .root_io import write_spectrum, read_spectrum

__all__ = [
    "ExpSpectrum",
    "EnergyCalibrator",
    "CalibrationResult",
    "GaussianFitter",
    "PeakFitter",
    "DeadTimeCorrection",
    "PileUpCorrection",
    "BackgroundSubtractor",
    "LivetimeDetermination",
    "LivetimeResult",
    "SpectrumPipeline",
    "PipelineResult",
    "process_spectrum",
    "process_run_sequence",
    "write_spectrum",
    "read_spectrum",
]
