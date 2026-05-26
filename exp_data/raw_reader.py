"""
Raw binary spectrum reader for 4096-channel MCA data.

Reads binary files with 128-byte header + 4096 channels of int32 counts.
Used for both Tc-99 beta spectra and Am-241 calibration sources.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.signal import find_peaks

from exp_data.spectrum import ExpSpectrum

# Binary format constants
HEADER_SIZE = 128  # bytes
CHANNEL_SIZE = 4  # int32
DEFAULT_N_CHANNELS = 4096
# 16512 = 128 + 4096 * 4
FILE_SIZE = HEADER_SIZE + DEFAULT_N_CHANNELS * CHANNEL_SIZE


def read_raw_spectrum(
    filepath: str | Path,
    n_channels: int = DEFAULT_N_CHANNELS,
    energy_keV: Optional[np.ndarray] = None,
    run_id: Optional[str] = None,
    source: Optional[str] = None,
    live_time_sec: Optional[float] = None,
    metadata: Optional[dict] = None,
) -> ExpSpectrum:
    """Read a raw binary spectrum file.

    Parameters
    ----------
    filepath
        Path to the binary spectrum file.
    n_channels
        Number of energy channels (default 4096).
    energy_keV
        Energy calibration array (n_channels,). If None, uses channel indices.
    run_id
        Run identifier string.
    source
        Source identifier (e.g. 'Tc99', 'Am241').
    live_time_sec
        Live acquisition time in seconds.
    metadata
        Additional metadata dict.

    Returns
    -------
    ExpSpectrum
        Parsed spectrum with metadata.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Spectrum file not found: {filepath}")

    file_size = filepath.stat().st_size
    expected = HEADER_SIZE + n_channels * CHANNEL_SIZE
    if file_size != expected:
        raise ValueError(
            f"File size {file_size} != expected {expected} "
            f"(header={HEADER_SIZE} + {n_channels}ch * {CHANNEL_SIZE}B)"
        )

    with open(filepath, "rb") as f:
        # Skip header
        f.seek(HEADER_SIZE)
        counts = np.fromfile(f, dtype=np.int32, count=n_channels)

    if energy_keV is None:
        energy_keV = np.arange(n_channels, dtype=np.float64)

    return ExpSpectrum(
        energies=energy_keV,
        counts=counts.astype(np.float64),
        metadata=metadata or {},
        run_id=run_id,
        source=source,
        live_time=live_time_sec,
    )


def subtract_spectra(
    cumulative_a: ExpSpectrum,
    cumulative_b: ExpSpectrum,
    run_id: Optional[str] = None,
) -> ExpSpectrum:
    """Subtract two cumulative spectra to get the incremental run.

    Both spectra must have the same number of channels.
    Result = A - B (where A is the later cumulative measurement).

    Parameters
    ----------
    cumulative_a
        Later cumulative spectrum (larger time).
    cumulative_b
        Earlier cumulative spectrum (smaller time).
    run_id
        Run ID for the difference spectrum.

    Returns
    -------
    ExpSpectrum
        The difference spectrum (A - B).
    """
    assert len(cumulative_a.counts) == len(cumulative_b.counts), (
        f"Spectra must have same channel count: "
        f"{len(cumulative_a.counts)} vs {len(cumulative_b.counts)}"
    )

    diff_counts = cumulative_a.counts - cumulative_b.counts
    # Propagate errors: σ_diff = sqrt(σ_a² + σ_b²)
    diff_errors = np.sqrt(cumulative_a.errors**2 + cumulative_b.errors**2)

    return ExpSpectrum(
        energies=cumulative_a.energies.copy(),
        counts=diff_counts,
        errors=diff_errors,
        metadata={
            **cumulative_a.metadata,
            "operation": "subtract",
            "parent_a": cumulative_a.run_id,
            "parent_b": cumulative_b.run_id,
        },
        run_id=run_id,
        source=cumulative_a.source,
    )


def calibrate_from_am241(
    filepath: str | Path,
    known_peaks: dict[float, int] | None = None,
    n_channels: int = DEFAULT_N_CHANNELS,
) -> tuple[float, float]:
    """Calibrate energy scale using Am-241 calibration source.

    Matches known Am-241 gamma lines to detected peaks in the
    low-energy region (ch < 300) and fits a linear calibration.

    Parameters
    ----------
    filepath
        Path to the Am-241 binary file (CAM1.DAT or CAM2.DAT).
    known_peaks
        Mapping of energy_keV -> expected channel. If None, auto-detects
        Am-241 peaks at 13.9, 17.8, 26.3 keV.
    n_channels
        Number of channels.

    Returns
    -------
    (offset_keV, gain_keV_per_ch)
        Linear calibration: E = offset + gain * channel
    """
    filepath = Path(filepath)
    with open(filepath, "rb") as f:
        f.read(128)  # skip header
        counts = np.array(
            struct.unpack(f"{n_channels}i", f.read(n_channels * 4)), dtype=np.int64
        )

    # Search for peaks in low-energy region only
    counts_low = counts.copy()
    counts_low[300:] = 0
    peaks, props = find_peaks(counts_low, height=10, distance=10)
    peak_energies = np.array(peaks, dtype=float)
    _peak_heights = props["peak_heights"]  # noqa: F841

    # Known Am-241 gamma lines (keV)
    am241_lines = {13.9: None, 17.8: None, 26.3: None}
    if known_peaks:
        am241_lines = known_peaks
    else:
        # Match lines to peaks by proximity
        for E in am241_lines:
            best = peak_energies[np.argmin(np.abs(peak_energies - E / 0.1))]
            am241_lines[E] = int(round(best))

    # Use all matched peaks for calibration
    energies = np.array(list(am241_lines.keys()))
    channels = np.array(list(am241_lines.values()))

    # Linear fit
    coeffs = np.polyfit(channels, energies, 1)
    gain = coeffs[0]
    offset = coeffs[1]

    return float(offset), float(gain)


def read_am241_calibration(
    filepath: str | Path,
    channel_energy_map: dict[int, float] | None = None,
    n_channels: int = DEFAULT_N_CHANNELS,
) -> ExpSpectrum:
    """Read Am-241 calibration spectrum.

    Parameters
    ----------
    filepath
        Path to the binary calibration file (CAM1.DAT or CAM2.DAT).
    channel_energy_map
        Mapping of channel index -> energy_keV for known peaks.
        If None, assumes linear calibration from channel indices.
    n_channels
        Number of channels.

    Returns
    -------
    ExpSpectrum
        Calibration spectrum.
    """
    filepath = Path(filepath)

    # Am-241 known peaks (keV): 13.9, 17.8, 26.3, 59.54
    # The spectrum typically has a few prominent peaks
    with open(filepath, "rb") as f:
        _header = f.read(128)
        counts_raw = struct.unpack(f"{n_channels}i", f.read(n_channels * 4))

    counts = np.array(counts_raw, dtype=np.float64)
    energy_keV = np.arange(n_channels, dtype=np.float64)

    metadata = {
        "source": "Am241",
        "calibration_type": "energy",
        "known_peaks_keV": [13.9, 17.8, 26.3, 59.54],
    }

    return ExpSpectrum(
        energies=energy_keV,
        counts=counts,
        metadata=metadata,
        source="Am241",
    )


def get_run_files(
    directory: str | Path,
    prefix: str = "A",
    sort_key: str = "natural",
) -> list[Path]:
    """Get sorted list of run files in directory.

    Parameters
    ----------
    directory
        Directory containing run files.
    prefix
        Filename prefix (default 'A').
    sort_key
        Sort method: 'natural' (A1, A2, ..., A10, A11) or 'lexicographic' (A1, A10, A100).

    Returns
    -------
    list[Path]
        Sorted file paths.
    """
    directory = Path(directory)
    runs = [f for f in directory.iterdir() if f.is_file() and f.name.startswith(prefix)]

    if sort_key == "natural":
        runs.sort(key=lambda p: _natural_sort_key(p.name))
    else:
        runs.sort(key=lambda p: p.name)

    return runs


def _natural_sort_key(s: str) -> list:
    """Sort key for natural ordering: A1, A2, ..., A10, A100."""
    import re

    parts = re.split(r"(\d+)", s)
    return [int(p) if p.isdigit() else p.lower() for p in parts]
