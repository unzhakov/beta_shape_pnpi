"""
Nuclear data integration module.

Provides interfaces to retrieve nuclear decay data from paceENSDF
and from custom JSON input files. Both sources produce a valid
SpectrumConfig for beta spectrum calculation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from beta_spectrum.spectrum import SpectrumConfig

# Mapping from paceENSDF forbiddenness codes to our transition_type
FORBIDDENNESS_MAP: Dict[str, str] = {
    "0A": "A",
    "1A": "A",
    "1F": "F1",
    "1UF": "F1U",
    "2F": "F2",
    "2UF": "F2U",
    "3F": "F3",
    "3UF": "F3U",
    "4F": "F4",
    "4UF": "F4",
    "5F": "F4",
    "5UF": "F4",
}


@dataclass
class DecayInfo:
    """
    Parsed nuclear decay information from paceENSDF or JSON.

    Attributes
    ----------
    parent_symbol : str
        Parent isotope symbol (e.g. "Tc99").
    Z_parent : int
        Atomic number of parent.
    Z_daughter : int
        Atomic number of daughter.
    A_number : int
        Mass number (same for parent and daughter in beta decay).
    endpoint_keV : float
        Q-value in keV.
    endpoint_MeV : float
        Q-value in MeV.
    transition_type : str
        Forbiddenness classification (A, F1, F1U, F2, F2U, ...).
    forbiddenness_code : str
        Original paceENSDF forbiddenness code (e.g. "2F").
    parent_J : float
        Parent spin.
    parent_pi : int
        Parent parity (1=+, -1=-).
    half_life : str
        Half-life string from ENSDF (e.g. "y", "h").
    branches : List[BranchInfo]
        List of decay branches to daughter states.
    """

    parent_symbol: str
    Z_parent: int
    Z_daughter: int
    A_number: int
    endpoint_keV: float
    endpoint_MeV: float
    transition_type: str
    forbiddenness_code: str
    parent_J: float
    parent_pi: int
    half_life: str
    branches: List[BranchInfo]


@dataclass
class BranchInfo:
    """
    A single decay branch to a daughter state.

    Attributes
    ----------
    level_index : int
        Level index in daughter nucleus.
    level_energy_keV : float
        Excitation energy of the daughter level.
    intensity : float
        Branch intensity (fraction of decays).
    log_ft : float
        log(ft) value for this branch.
    """

    level_index: int
    level_energy_keV: float
    intensity: float
    log_ft: Optional[float]


def _parse_nuclide_symbol(symbol: str) -> Tuple[str, int, int]:
    """
    Parse a nuclide symbol like 'Tc99' into (element, Z, A).

    Parameters
    ----------
    symbol : str
        Nuclide symbol, e.g. "Tc99", "Co60", "Y86".

    Returns
    -------
    tuple[str, int, int]
        (element_symbol, Z, A)

    Raises
    ------
    ValueError
        If the symbol cannot be parsed.
    """
    # Separate letters from digits
    i = 0
    while i < len(symbol) and symbol[i].isalpha():
        i += 1

    if i == 0 or i == len(symbol):
        raise ValueError(
            f"Cannot parse nuclide symbol '{symbol}'. "
            f"Expected format like 'Tc99', 'Co60', 'Y86'."
        )

    element = symbol[:i]
    A = int(symbol[i:])

    # Look up Z from element symbol
    element_to_Z: Dict[str, int] = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
        "Cs": 55,
        "Ba": 56,
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "At": 85,
        "Rn": 86,
        "Fr": 87,
        "Ra": 88,
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        "Np": 93,
        "Pu": 94,
        "Am": 95,
        "Cm": 96,
        "Bk": 97,
        "Cf": 98,
    }

    if element not in element_to_Z:
        raise ValueError(
            f"Unknown element symbol '{element}'. "
            f"Supported: {sorted(element_to_Z.keys())}"
        )

    return element, element_to_Z[element], A


def _get_decay_mode_symbol(
    decay_type: str,
) -> str:
    """
    Map internal decay type to paceENSDF mode string.

    Parameters
    ----------
    decay_type : str
        One of 'beta_minus', 'beta_plus', 'ec'.

    Returns
    -------
    str
        paceENSDF mode string: 'BM', 'ECBP'.
    """
    mapping = {
        "beta_minus": "BM",
        "beta_plus": "ECBP",
        "ec": "ECBP",
    }
    if decay_type not in mapping:
        raise ValueError(
            f"Unknown decay type '{decay_type}'. " f"Supported: {list(mapping.keys())}"
        )
    return mapping[decay_type]


def _resolve_decay_index(
    edata: list[Any],
    parent_symbol: str,
    mode: str,
    decay_index: Optional[int] = None,
) -> int:
    """
    Find the correct decay index for a parent isotope.

    Parameters
    ----------
    edata : list
        paceENSDF data from load_ensdf().
    parent_symbol : str
        Parent isotope symbol.
    mode : str
        paceENSDF mode ('BM', 'ECBP', 'A').
    decay_index : int, optional
        Explicit index. If None, selects ground state (0) if available,
        otherwise the first available decay.

    Returns
    -------
    int
        The decay index to use.

    Raises
    ------
    ValueError
        If no matching decay is found.
    """
    import paceENSDF as pe

    e = pe.ENSDF()
    pairs = e.ensdf_pairs(edata, mode)

    matching = [(k, v) for k, v in pairs.items() if str(k[0]) == parent_symbol]

    if not matching:
        raise ValueError(
            f"No {mode} decay data found for {parent_symbol}. "
            f"Available {mode} parents: {[str(k[0]) for k in pairs.keys()]}"
        )

    if decay_index is not None:
        for k, v in matching:
            if k[3] == decay_index:
                return decay_index
        raise ValueError(
            f"Decay index {decay_index} not found for {parent_symbol}. "
            f"Available indices: {[k[3] for k, v in matching]}"
        )

    # Default: ground state (index 0), or first available
    for k, v in matching:
        if k[3] == 0:
            return 0
    return int(matching[0][0][3])


def get_decay_info_from_paceENSDF(
    nuclide: str,
    decay_type: str = "beta_minus",
    decay_index: Optional[int] = None,
) -> DecayInfo:
    """
    Retrieve decay information for a nuclide from paceENSDF.

    Parameters
    ----------
    nuclide : str
        Nuclide symbol, e.g. "Tc99", "Co60".
    decay_type : str
        Decay mode: 'beta_minus', 'beta_plus', or 'ec'.
    decay_index : int, optional
        Decay index (0 = ground state). If None, defaults to ground state.

    Returns
    -------
    DecayInfo
        Parsed decay information.

    Examples
    --------
    >>> info = get_decay_info_from_paceENSDF("Tc99", "beta_minus")
    >>> config = decay_info_to_config(info)
    """
    import paceENSDF as pe

    e = pe.ENSDF()
    edata = e.load_ensdf()

    mode = _get_decay_mode_symbol(decay_type)
    _element, Z_parent, A_number = _parse_nuclide_symbol(nuclide)
    decay_idx = _resolve_decay_index(edata, nuclide, mode, decay_index)

    # Get parent decay info
    _decay_data = e.get_parent_decay(edata, nuclide, decay_idx, mode=mode)
    halflife_data = e.get_parent_halflife(
        edata, nuclide, decay_idx, mode=mode, units="best"
    )
    jpi_data = e.get_parent_jpi(edata, nuclide, decay_idx, mode=mode)
    beta_data = e.get_beta_minus(edata, nuclide, decay_idx, units="best")

    # Extract Q-value and half-life
    q_value_keV = 0.0
    half_life_str = "unknown"
    for k, v in halflife_data.items():
        q_value_keV = v[0]
        half_life_str = str(v[2]) if len(v) > 2 else "unknown"
        break

    # Extract parent spin/parity
    parent_J = jpi_data[0][0] if jpi_data else 0.0
    parent_pi = jpi_data[0][1] if jpi_data else 1

    # Extract daughter Z from decay data
    Z_daughter = Z_parent + 1 if decay_type == "beta_minus" else Z_parent - 1

    # Extract branches and forbiddenness info
    branches: List[BranchInfo] = []
    forbidden_code = "0A"  # default
    transition_type = "A"  # default

    for k, v in beta_data.items():
        for state in v:
            # state[15] contains forbiddenness code (e.g. '2F', '1UF')
            if len(state) > 15 and state[15] is not None:
                code = str(state[15])
                if code in FORBIDDENNESS_MAP:
                    forbidden_code = code
                    transition_type = FORBIDDENNESS_MAP[code]

            branch = BranchInfo(
                level_index=int(state[0]),
                level_energy_keV=float(state[1]),
                intensity=float(state[8]),
                log_ft=float(state[11]) if state[11] is not None else None,
            )
            branches.append(branch)

    return DecayInfo(
        parent_symbol=nuclide,
        Z_parent=Z_parent,
        Z_daughter=Z_daughter,
        A_number=A_number,
        endpoint_keV=q_value_keV,
        endpoint_MeV=q_value_keV / 1000.0,
        transition_type=transition_type,
        forbiddenness_code=forbidden_code,
        parent_J=parent_J,
        parent_pi=parent_pi,
        half_life=half_life_str,
        branches=branches,
    )


def decay_info_to_config(
    info: DecayInfo,
    e_step_MeV: float = 0.001,
    use_detector_response: bool = False,
) -> SpectrumConfig:
    """
    Convert DecayInfo to a SpectrumConfig.

    Parameters
    ----------
    info : DecayInfo
        Parsed decay information.
    e_step_MeV : float
        Energy step size in MeV.
    use_detector_response : bool
        Whether to enable detector response convolution.

    Returns
    -------
    SpectrumConfig
        Configuration ready for BetaSpectrum.from_config().
    """
    return SpectrumConfig(
        Z_parent=info.Z_parent,
        Z_daughter=info.Z_daughter,
        A_number=info.A_number,
        endpoint_MeV=info.endpoint_MeV,
        transition_type=info.transition_type,
        e_step_MeV=e_step_MeV,
        use_detector_response=use_detector_response,
    )


# ---------------------------------------------------------------------------
# JSON input support
# ---------------------------------------------------------------------------

# Default JSON schema
DEFAULT_JSON_SCHEMA = {
    "parent_symbol": "Tc99",
    "decay_type": "beta_minus",
    "endpoint_MeV": 0.2111,
    "Z_parent": 43,
    "Z_daughter": 44,
    "A_number": 99,
    "transition_type": "F2",
    "e_step_MeV": 0.001,
    "use_detector_response": False,
    "detector_model": "gaussian",
    "detector_sigma_a_keV": 1.0,
    "detector_sigma_b": 0.0,
    "detector_tail_fraction": 0.0,
    "detector_tau_keV": 5.0,
    "detector_fano_factor": 0.12,
    "detector_n_channels": 4096,
    "detector_channel_energy_range": [0.0, 0.35],
}


def validate_json_input(data: Dict[str, Any]) -> None:
    """
    Validate a JSON input dictionary.

    Parameters
    ----------
    data : dict
        Parsed JSON data.

    Raises
    ------
    ValueError
        If required fields are missing or invalid.
    """
    required = ["endpoint_MeV", "Z_parent", "Z_daughter", "A_number"]
    for field in required:
        if field not in data:
            raise ValueError(f"Missing required field: '{field}'")

    if data["endpoint_MeV"] <= 0:
        raise ValueError("endpoint_MeV must be positive")

    if data["Z_parent"] <= 0 or data["Z_daughter"] <= 0:
        raise ValueError("Z_parent and Z_daughter must be positive")

    if data["A_number"] <= 0:
        raise ValueError("A_number must be positive")

    if data.get("transition_type", "A") not in [
        "A",
        "F1",
        "F1U",
        "F2",
        "F2U",
        "F3",
        "F3U",
        "F4",
    ]:
        raise ValueError(
            f"Invalid transition_type: '{data['transition_type']}'. "
            f"Supported: A, F1, F1U, F2, F2U, F3, F3U, F4"
        )


def load_json_input(filepath: str) -> Dict[str, Any]:
    """
    Load and validate a JSON input file.

    Parameters
    ----------
    filepath : str
        Path to the JSON file.

    Returns
    -------
    dict
        Validated input data with defaults filled in.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    with open(path, "r") as f:
        data = json.load(f)

    # Fill in defaults for missing optional fields
    for key, default in DEFAULT_JSON_SCHEMA.items():
        if key not in data:
            data[key] = default

    validate_json_input(data)
    return data


def json_to_config(data: Dict[str, Any]) -> SpectrumConfig:
    """
    Convert a validated JSON input dictionary to SpectrumConfig.

    Parameters
    ----------
    data : dict
        Validated JSON input data.

    Returns
    -------
    SpectrumConfig
        Configuration ready for BetaSpectrum.from_config().
    """
    return SpectrumConfig(
        Z_parent=int(data["Z_parent"]),
        Z_daughter=int(data["Z_daughter"]),
        A_number=int(data["A_number"]),
        endpoint_MeV=float(data["endpoint_MeV"]),
        transition_type=str(data.get("transition_type", "A")),
        e_step_MeV=float(data.get("e_step_MeV", 0.001)),
        use_detector_response=bool(data.get("use_detector_response", False)),
        detector_model=str(data.get("detector_model", "gaussian")),
        detector_sigma_a_keV=float(data.get("detector_sigma_a_keV", 1.0)),
        detector_sigma_b=float(data.get("detector_sigma_b", 0.0)),
        detector_tail_fraction=float(data.get("detector_tail_fraction", 0.0)),
        detector_tau_keV=float(data.get("detector_tau_keV", 5.0)),
        detector_fano_factor=float(data.get("detector_fano_factor", 0.12)),
        detector_n_channels=int(data.get("detector_n_channels", 4096)),
        detector_channel_energy_range=tuple(
            data.get("detector_channel_energy_range", [0.0, 0.35])
        ),
    )


def create_config_from_source(
    source: str,
    nuclide: Optional[str] = None,
    json_path: Optional[str] = None,
    **kwargs: Any,
) -> SpectrumConfig:  # type: ignore[return-value]
    """
    Create a SpectrumConfig from either paceENSDF or JSON input.

    Parameters
    ----------
    source : str
        Data source: 'paceENSDF' or 'json'.
    nuclide : str, optional
        Nuclide symbol for paceENSDF source (e.g. "Tc99").
    json_path : str, optional
        Path to JSON input file for JSON source.
    **kwargs
        Additional parameters passed to decay_info_to_config or json_to_config.

    Returns
    -------
    SpectrumConfig
        Configuration ready for BetaSpectrum.from_config().

    Raises
    ------
    ValueError
        If source is unknown or required arguments are missing.

    Examples
    --------
    >>> config = create_config_from_source("paceENSDF", nuclide="Tc99")
    >>> config = create_config_from_source("json", json_path="custom_input.json")
    """
    if source == "paceENSDF":
        if not nuclide:
            raise ValueError("nuclide parameter required for paceENSDF source")
        decay_type = kwargs.get("decay_type", "beta_minus")
        decay_index = kwargs.get("decay_index", None)
        e_step = kwargs.get("e_step_MeV", 0.001)
        use_det = kwargs.get("use_detector_response", False)
        info = get_decay_info_from_paceENSDF(nuclide, decay_type, decay_index)
        return decay_info_to_config(
            info, e_step_MeV=e_step, use_detector_response=use_det
        )
    elif source == "json":
        if not json_path:
            raise ValueError("json_path parameter required for JSON source")
        data = load_json_input(json_path)
        return json_to_config(data)
    else:
        raise ValueError(f"Unknown source: '{source}'. Supported: 'paceENSDF', 'json'")
