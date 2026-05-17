"""
Fitter module — main analysis framework for beta spectrometry.

Provides the analysis workhorse: theoretical model + detector convolution,
multi-parameter fitting, and physics parameter extraction (C(W), g_V, g_A).

Usage:
    from fitter import SpectrumModel, SpectrumFitter, AnalysisResult

    # Build theoretical model with detector response
    model = SpectrumModel(spectrum, detector_response, channel_energies)

    # Fit model to experimental data
    fitter = SpectrumFitter(model, exp_spectrum)
    result = fitter.fit(
        x0=[endpoint_keV, normalization, background],
        bounds=[(...), (...), (...)],
    )

    # Extract physics parameters
    print(result.summary())
"""

from .model import SpectrumModel
from .fit_engine import SpectrumFitter, CurveFitter, FitConfig
from .result import AnalysisResult, FitSummary
from .extractor import CWExtractor, GVAExtractor

__all__ = [
    "SpectrumModel",
    "SpectrumFitter",
    "CurveFitter",
    "FitConfig",
    "AnalysisResult",
    "FitSummary",
    "CWExtractor",
    "GVAExtractor",
]
