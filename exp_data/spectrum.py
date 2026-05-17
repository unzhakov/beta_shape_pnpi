"""
Experimental spectrum container.

Stores measured detector spectra with metadata, uncertainties,
and calibration information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class ExpSpectrum:
    """Container for an experimental beta spectrum measurement.

    Parameters
    ----------
    energies : np.ndarray
        Energy calibration values for each channel (keV), shape (N_channels,).
    counts : np.ndarray
        Raw counts per channel, shape (N_channels,).
    errors : np.ndarray, optional
        Statistical uncertainties per channel. Defaults to sqrt(counts).
    metadata : dict, optional
        Arbitrary metadata (source, run ID, date, etc.).
    dead_time : float, optional
        Dead time fraction (0.0 to 1.0).
    live_time : float, optional
        Live acquisition time in seconds.
    source : str, optional
        Radioactive source identifier.
    run_id : str, optional
        Unique run identifier.
    date : str, optional
        Measurement date (ISO format).
    """

    energies: np.ndarray
    counts: np.ndarray
    errors: Optional[np.ndarray] = None
    metadata: dict = field(default_factory=dict)
    dead_time: float = 0.0
    live_time: Optional[float] = None
    source: Optional[str] = None
    run_id: Optional[str] = None
    date: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate and compute default errors."""
        self.energies = np.asarray(self.energies, dtype=np.float64)
        self.counts = np.asarray(self.counts, dtype=np.float64)
        assert self.energies.shape == self.counts.shape, (
            f"Energies ({self.energies.shape}) and counts ({self.counts.shape}) "
            "must have the same shape"
        )
        assert np.all(self.energies >= 0), "Energies must be non-negative"
        assert np.all(self.counts >= 0), "Counts must be non-negative"

        if self.errors is None:
            # Poisson statistics: σ = √N (handle zero counts gracefully)
            self.errors = np.sqrt(np.maximum(self.counts, 1.0))
        else:
            self.errors = np.asarray(self.errors, dtype=np.float64)
            assert self.errors.shape == self.counts.shape
            assert np.all(self.errors > 0), "Errors must be positive"

    @property
    def n_channels(self) -> int:
        """Number of energy channels."""
        return len(self.energies)

    @property
    def energy_range(self) -> tuple[float, float]:
        """(min_energy_keV, max_energy_keV)."""
        return (float(self.energies[0]), float(self.energies[-1]))

    @property
    def total_counts(self) -> float:
        """Total counts in the spectrum."""
        return float(np.sum(self.counts))

    @property
    def live_time_corrected_rate(self) -> Optional[float]:
        """Count rate corrected for dead time, in counts per second."""
        if self.live_time is None or self.live_time <= 0:
            return None
        live_seconds = self.live_time * (1.0 - self.dead_time)
        return self.total_counts / live_seconds

    def normalize(self) -> np.ndarray:
        """Return spectrum normalized to unit area."""
        integral = np.trapezoid(self.counts, self.energies)
        if integral > 0:
            return self.counts / integral
        return self.counts.copy()

    def apply_dead_time_correction(self) -> ExpSpectrum:
        """Return a new spectrum corrected for dead time.

        Correction: N_true = N_measured / (1 - τ), where τ is dead time fraction.
        """
        if self.dead_time <= 0 or self.dead_time >= 1.0:
            return ExpSpectrum(
                energies=self.energies.copy(),
                counts=self.counts.copy(),
                errors=self.errors.copy(),
                metadata=self.metadata.copy(),
                dead_time=self.dead_time,
                live_time=self.live_time,
                source=self.source,
                run_id=self.run_id,
                date=self.date,
            )

        corrected_counts = self.counts / (1.0 - self.dead_time)
        corrected_errors = self.errors / (1.0 - self.dead_time)

        return ExpSpectrum(
            energies=self.energies.copy(),
            counts=corrected_counts,
            errors=corrected_errors,
            metadata=self.metadata.copy(),
            dead_time=0.0,  # Corrected
            live_time=self.live_time,
            source=self.source,
            run_id=self.run_id,
            date=self.date,
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "energies": self.energies.tolist(),
            "counts": self.counts.tolist(),
            "errors": self.errors.tolist(),
            "metadata": self.metadata,
            "dead_time": self.dead_time,
            "live_time": self.live_time,
            "source": self.source,
            "run_id": self.run_id,
            "date": self.date,
        }

    @classmethod
    def from_dict(cls, data: dict) -> ExpSpectrum:
        """Deserialize from dictionary."""
        return cls(
            energies=np.array(data["energies"], dtype=np.float64),
            counts=np.array(data["counts"], dtype=np.float64),
            errors=np.array(data["errors"], dtype=np.float64),
            metadata=data.get("metadata", {}),
            dead_time=data.get("dead_time", 0.0),
            live_time=data.get("live_time"),
            source=data.get("source"),
            run_id=data.get("run_id"),
            date=data.get("date"),
        )
