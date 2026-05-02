# beta_spectrum/__init__.py

from .constants import ALPHA, HBAR_C_MEV_FM, ME_MEV, MP_MEV
from .utils import T_to_W, W_to_T, momentum

from .components.phase_space import PhaseSpace
from .components.fermi import FermiFunction
from .components.finite_size import FiniteSizeL0, ChargeDistributionU
from .components.screening import ScreeningCorrection
from .components.exchange import ExchangeCorrection
from .components.radiative import RadiativeCorrection
from .components.detector_response import DetectorResponse

from .spectrum import BetaSpectrum, SpectrumConfig, BetaSpectrumAnalyzer
from .fitter import CurveFitter, FitConfig, FitResult
from .cw_extractor import CWExtractor, CWExtractionResult, GVAExtractionResult
from .nuclear_data import (
    DecayInfo,
    BranchInfo,
    get_decay_info_from_paceENSDF,
    decay_info_to_config,
    create_config_from_source,
    load_json_input,
    json_to_config,
    DEFAULT_JSON_SCHEMA,
)

__all__ = [
    "ALPHA",
    "HBAR_C_MEV_FM",
    "ME_MEV",
    "MP_MEV",
    "T_to_W",
    "W_to_T",
    "momentum",
    "PhaseSpace",
    "FermiFunction",
    "FiniteSizeL0",
    "ChargeDistributionU",
    "ScreeningCorrection",
    "ExchangeCorrection",
    "RadiativeCorrection",
    "DetectorResponse",
    "BetaSpectrum",
    "SpectrumConfig",
    "BetaSpectrumAnalyzer",
    "CurveFitter",
    "FitConfig",
    "FitResult",
    "CWExtractor",
    "CWExtractionResult",
    "GVAExtractionResult",
    "DecayInfo",
    "BranchInfo",
    "get_decay_info_from_paceENSDF",
    "decay_info_to_config",
    "create_config_from_source",
    "load_json_input",
    "json_to_config",
    "DEFAULT_JSON_SCHEMA",
]

__version__ = "0.3.0"
