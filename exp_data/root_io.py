"""
ROOT file I/O for experimental spectra.

Supports PyROOT (primary) and uproot (fallback) for reading and writing
ExpSpectrum objects to ROOT TTree/TBranch format.

Usage:
    from exp_data.root_io import write_spectrum, read_spectrum

    # Write
    write_spectrum("output.root", spectrum)

    # Read
    spectrum = read_spectrum("output.root")
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from exp_data.spectrum import ExpSpectrum


# Try PyROOT first, fall back to uproot
_PYROOT_AVAILABLE = False
_UPTOOL_AVAILABLE = False

try:
    import ROOT  # noqa: F401

    _PYROOT_AVAILABLE = True
except ImportError:
    try:
        import uproot  # type: ignore[import-not-found]  # noqa: F401

        _UPTOOL_AVAILABLE = True
    except ImportError:
        raise ImportError(
            "Neither PyROOT nor uproot is available. "
            "Install one of: pip install beta-spectrum[root] (PyROOT) "
            "or pip install uproot (pure Python fallback)."
        )


def write_spectrum(
    filepath: str | Path,
    spectrum: "ExpSpectrum",
    tree_name: str = "spectrum",
) -> None:
    """Write an ExpSpectrum to a ROOT file.

    Creates a TTree with TBranches for energies, counts, errors,
    scalar fields (dead_time, live_time, source, run_id, date),
    and metadata as a JSON string.

    Parameters
    ----------
    filepath : str or Path
        Output ROOT file path.
    spectrum : ExpSpectrum
        Spectrum to write.
    tree_name : str
        Name of the TTree in the ROOT file.
    """
    filepath = Path(filepath)

    if _PYROOT_AVAILABLE:
        _write_pyroot(filepath, spectrum, tree_name)
    elif _UPTOOL_AVAILABLE:
        _write_uproot(filepath, spectrum, tree_name)
    else:
        raise RuntimeError("No ROOT I/O library available")


def read_spectrum(
    filepath: str | Path,
    tree_name: str = "spectrum",
) -> "ExpSpectrum":
    """Read an ExpSpectrum from a ROOT file.

    Parameters
    ----------
    filepath : str or Path
        Input ROOT file path.
    tree_name : str
        Name of the TTree in the ROOT file.

    Returns
    -------
    ExpSpectrum
        Reconstructed spectrum with all metadata.
    """
    filepath = Path(filepath)

    if _PYROOT_AVAILABLE:
        return _read_pyroot(filepath, tree_name)
    elif _UPTOOL_AVAILABLE:
        return _read_uproot(filepath, tree_name)
    else:
        raise RuntimeError("No ROOT I/O library available")


# --- PyROOT implementation ---


def _write_pyroot(
    filepath: Path,
    spectrum: "ExpSpectrum",
    tree_name: str,
) -> None:
    """Write ExpSpectrum using PyROOT TTree/TBranch."""
    import ROOT
    from array import array

    with ROOT.TFile(str(filepath), "RECREATE"):
        tree = ROOT.TTree(tree_name, "Experimental beta spectrum")

        # Data arrays
        energies_arr = np.array(spectrum.energies, dtype=np.float64)
        counts_arr = np.array(spectrum.counts, dtype=np.float64)
        errors_arr = np.array(spectrum.errors, dtype=np.float64)

        # TBranches for arrays
        n = len(energies_arr)
        tree.Branch("energies", energies_arr, f"energies[{n}]/D")
        tree.Branch("counts", counts_arr, f"counts[{n}]/D")
        tree.Branch("errors", errors_arr, f"errors[{n}]/D")

        # Scalar TBranches (using array of length 1)
        dead_time_val = array("d", [spectrum.dead_time])
        live_time_val = array(
            "d", [spectrum.live_time if spectrum.live_time is not None else 0.0]
        )

        tree.Branch("dead_time", dead_time_val, "dead_time/D")
        tree.Branch("live_time", live_time_val, "live_time/D")

        # String fields
        source_str = ROOT.TObjString(spectrum.source or "")
        run_id_str = ROOT.TObjString(spectrum.run_id or "")
        date_str = ROOT.TObjString(spectrum.date or "")

        tree.Branch("source", source_str)
        tree.Branch("run_id", run_id_str)
        tree.Branch("date", date_str)

        # Metadata as JSON
        metadata_json = json.dumps(spectrum.metadata)
        metadata_str = ROOT.TObjString(metadata_json)
        tree.Branch("metadata_json", metadata_str)

        tree.Fill()
        tree.Write()


def _read_pyroot(
    filepath: Path,
    tree_name: str,
) -> "ExpSpectrum":
    """Read ExpSpectrum from ROOT file using PyROOT."""
    import ROOT
    import json

    f = ROOT.TFile(str(filepath), "READ")
    tree = f.Get(tree_name)

    # Read first (and only) entry
    tree.GetEntry(0)

    energies = np.array(tree.energies, dtype=np.float64)
    counts = np.array(tree.counts, dtype=np.float64)
    errors = np.array(tree.errors, dtype=np.float64)

    dead_time = float(tree.dead_time)
    live_time = float(tree.live_time) if tree.live_time > 0 else None
    source_raw = tree.source
    source = str(source_raw) if source_raw and str(source_raw) else None
    run_id_raw = tree.run_id
    run_id = str(run_id_raw) if run_id_raw and str(run_id_raw) else None
    date_raw = tree.date
    date = str(date_raw) if date_raw and str(date_raw) else None

    # Parse metadata JSON
    metadata_raw = tree.metadata_json
    # Handle various ROOT string types (TObjString, TString, str, bytes)
    if hasattr(metadata_raw, "__str__"):
        try:
            metadata_str = str(metadata_raw)
        except Exception:
            metadata_str = ""
    elif hasattr(metadata_raw, "GetString"):
        metadata_str = str(metadata_raw.GetString())
    else:
        metadata_str = ""
    metadata = json.loads(metadata_str) if metadata_str else {}

    from exp_data.spectrum import ExpSpectrum

    return ExpSpectrum(
        energies=energies,
        counts=counts,
        errors=errors,
        metadata=metadata,
        dead_time=dead_time,
        live_time=live_time,
        source=source,
        run_id=run_id,
        date=date,
    )


# --- uproot implementation ---


def _write_uproot(
    filepath: Path,
    spectrum: "ExpSpectrum",
    tree_name: str,
) -> None:
    """Write ExpSpectrum using uproot (pure Python)."""
    import uproot

    with uproot.recreate(str(filepath)) as f:
        f[tree_name] = {
            "energies": np.array(spectrum.energies, dtype=">f8"),
            "counts": np.array(spectrum.counts, dtype=">f8"),
            "errors": np.array(spectrum.errors, dtype=">f8"),
            "dead_time": np.array([spectrum.dead_time], dtype=">f8"),
            "live_time": np.array(
                [spectrum.live_time if spectrum.live_time is not None else 0.0],
                dtype=">f8",
            ),
            "source": np.array([spectrum.source or ""], dtype="S100"),
            "run_id": np.array([spectrum.run_id or ""], dtype="S100"),
            "date": np.array([spectrum.date or ""], dtype="S100"),
            "metadata_json": np.array([json.dumps(spectrum.metadata)], dtype="S500"),
        }


def _read_uproot(
    filepath: Path,
    tree_name: str,
) -> "ExpSpectrum":
    """Read ExpSpectrum from ROOT file using uproot."""
    import uproot
    import json

    with uproot.open(str(filepath)) as f:
        tree = f[tree_name]
        data = tree.arrays(library="np")

        energies = np.array(data["energies"], dtype=np.float64)
        counts = np.array(data["counts"], dtype=np.float64)
        errors = np.array(data["errors"], dtype=np.float64)

        dead_time = float(data["dead_time"][0])
        live_time_raw = float(data["live_time"][0])
        live_time = live_time_raw if live_time_raw > 0 else None

        source_raw = data["source"][0]
        source = source_raw.decode() if isinstance(source_raw, bytes) else source_raw
        run_id_raw = data["run_id"][0]
        run_id = run_id_raw.decode() if isinstance(run_id_raw, bytes) else run_id_raw
        date_raw = data["date"][0]
        date = date_raw.decode() if isinstance(date_raw, bytes) else date_raw

        metadata_str = data["metadata_json"][0]
        if isinstance(metadata_str, bytes):
            metadata_str = metadata_str.decode()
        metadata = json.loads(metadata_str) if metadata_str else {}

        return ExpSpectrum(
            energies=energies,
            counts=counts,
            errors=errors,
            metadata=metadata,
            dead_time=dead_time,
            live_time=live_time,
            source=source if source else None,
            run_id=run_id if run_id else None,
            date=date if date else None,
        )
