"""
Experimental data processing pipeline.

Orchestrates the workflow: read -> subtract -> calibrate -> livetime -> validate -> export.

Usage:
    from exp_data.pipeline import process_spectrum, process_run_sequence

    # Single spectrum
    result = process_spectrum(
        raw_path="data/raw/A1.DAT",
        calibration_path="data/raw/CAM1.DAT",
        pulser_frequency_hz=100.0,
        real_time_sec=3600.0,
    )

    # Run sequence (cumulative subtraction)
    results = process_run_sequence(
        run_dir="data/raw/",
        prefix="A",
        calibration_path="data/raw/CAM1.DAT",
        pulser_frequency_hz=100.0,
        real_time_sec=3600.0,
    )
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from exp_data.spectrum import ExpSpectrum
from exp_data.raw_reader import read_raw_spectrum, subtract_spectra, get_run_files
from exp_data.calibration import EnergyCalibrator, CalibrationResult
from exp_data.livetime import LivetimeDetermination, LivetimeResult


@dataclass
class PipelineResult:
    """Result of a pipeline execution."""

    spectrum: ExpSpectrum
    """Calibrated, corrected ExpSpectrum."""

    calibration: CalibrationResult
    """Energy calibration result."""

    livetime: Optional[LivetimeResult]
    """Live-time determination result (None if not computed)."""

    validation_passed: bool
    """Whether all quality checks passed."""

    validation_messages: list[str] = field(default_factory=list)
    """Validation messages (warnings/errors)."""


class SpectrumPipeline:
    """Orchestrate the full experimental data processing pipeline.

    Workflow:
        1. Read raw binary spectrum
        2. (Optional) Subtract from previous cumulative spectrum
        3. Apply energy calibration from Am-241
        4. Determine live-time from pulser (if available)
        5. Validate quality
        6. Export (optional)
    """

    def __init__(
        self,
        calibration_path: str | Path,
        pulser_frequency_hz: Optional[float] = None,
        real_time_sec: Optional[float] = None,
    ):
        self.calibration_path = Path(calibration_path)
        self.pulser_frequency_hz = pulser_frequency_hz
        self.real_time_sec = real_time_sec

    def process(
        self,
        raw_path: str | Path,
        cumulative_path: Optional[str | Path] = None,
        run_id: Optional[str] = None,
        source: Optional[str] = None,
    ) -> PipelineResult:
        """Process a single spectrum through the full pipeline.

        Parameters
        ----------
        raw_path : str or Path
            Path to the raw binary spectrum file.
        cumulative_path : str or Path, optional
            Path to the previous cumulative spectrum for subtraction.
        run_id : str, optional
            Run identifier.
        source : str, optional
            Source identifier.

        Returns
        -------
        PipelineResult
        """
        # Step 1: Read raw spectrum
        spectrum = read_raw_spectrum(
            raw_path,
            run_id=run_id,
            source=source,
            metadata={"operation": "read"},
        )

        # Step 2: Subtract from cumulative (if provided)
        if cumulative_path is not None:
            cumulative = read_raw_spectrum(cumulative_path)
            spectrum = subtract_spectra(spectrum, cumulative, run_id=run_id)

        # Step 3: Apply energy calibration
        calibration = EnergyCalibrator.from_am241_spectrum(self.calibration_path)
        spectrum = spectrum.apply_calibration(calibration)

        # Step 4: Determine live-time (if pulser data available)
        livetime = None
        if self.pulser_frequency_hz is not None and self.real_time_sec is not None:
            livetime = LivetimeDetermination.from_file(
                self.calibration_path,
                pulser_frequency_hz=self.pulser_frequency_hz,
                real_time_sec=self.real_time_sec,
                run_id=f"{run_id}_livetime" if run_id else "livetime",
                source="pulser",
            )
            # Store live-time fraction in metadata
            spectrum.metadata["live_time"] = {
                "fraction": livetime.live_time_fraction,
                "real_time_sec": livetime.real_time_sec,
                "pulser_frequency_hz": livetime.pulser_frequency_hz,
                "observed_pulser_rate": livetime.observed_pulser_rate,
                "method": "pulser",
            }
            spectrum.live_time = livetime.real_time_sec * livetime.live_time_fraction

        # Step 5: Validate quality
        messages = self._validate(spectrum, calibration, livetime)
        passed = len([m for m in messages if m.startswith("ERROR:")]) == 0

        return PipelineResult(
            spectrum=spectrum,
            calibration=calibration,
            livetime=livetime,
            validation_passed=passed,
            validation_messages=messages,
        )

    def process_sequence(
        self,
        run_dir: str | Path,
        prefix: str = "A",
        first_file: Optional[str | Path] = None,
        run_id: Optional[str] = None,
        source: Optional[str] = None,
    ) -> list[PipelineResult]:
        """Process a sequence of cumulative binary files.

        Reads all files matching prefix in run_dir, subtracts consecutive
        pairs to extract per-run spectra.

        Parameters
        ----------
        run_dir : str or Path
            Directory containing run files.
        prefix : str
            Filename prefix (default 'A').
        first_file : str or Path, optional
            Path to the first cumulative file (baseline for subtraction).
            If None, the first file in the sorted list is used as baseline.
        run_id : str, optional
            Base run identifier (appended with run number).
        source : str, optional
            Source identifier.

        Returns
        -------
        list of PipelineResult
            One result per extracted run.
        """
        files = get_run_files(run_dir, prefix=prefix)

        if len(files) < 2:
            raise ValueError(
                f"Need at least 2 cumulative files for subtraction, "
                f"found {len(files)} in {run_dir}"
            )

        results = []
        for i, filepath in enumerate(files[1:], start=1):
            result = self.process(
                raw_path=filepath,
                cumulative_path=files[0],
                run_id=f"{run_id}_{i}" if run_id else None,
                source=source,
            )
            results.append(result)

        return results

    def _validate(
        self,
        spectrum: ExpSpectrum,
        calibration: CalibrationResult,
        livetime: Optional[LivetimeResult],
    ) -> list[str]:
        """Validate spectrum quality.

        Checks:
        - Calibration fit quality (chi2/dof < 3)
        - Minimum total counts (> 100)
        - Live-time fraction in valid range (if available)
        - Energy range sanity (energies increasing, positive)
        """
        messages = []

        # Calibration quality
        if calibration.chi2_per_dof > 3.0:
            messages.append(
                f"ERROR: Calibration fit quality poor (chi2/dof = {calibration.chi2_per_dof:.2f} > 3.0)"
            )
        elif calibration.chi2_per_dof > 1.5:
            messages.append(
                f"WARNING: Calibration fit quality marginal (chi2/dof = {calibration.chi2_per_dof:.2f})"
            )

        # Minimum counts
        total = spectrum.total_counts
        if total < 100:
            messages.append(
                f"ERROR: Insufficient statistics (total counts = {total:.0f} < 100)"
            )

        # Live-time fraction
        if livetime is not None:
            if livetime.live_time_fraction < 0.5:
                messages.append(
                    f"WARNING: Low live-time fraction (tau = {livetime.live_time_fraction:.3f} < 0.5)"
                )
            elif livetime.live_time_fraction <= 0.0:
                messages.append(
                    f"ERROR: Invalid live-time fraction (tau = {livetime.live_time_fraction:.3f} <= 0)"
                )

        # Energy range sanity
        if spectrum.energies[0] < 0:
            messages.append("ERROR: Negative minimum energy")
        if not np.all(np.diff(spectrum.energies) >= 0):
            messages.append("ERROR: Energies not monotonically increasing")

        return messages


# Convenience functions


def process_spectrum(
    raw_path: str | Path,
    calibration_path: str | Path,
    run_id: Optional[str] = None,
    source: Optional[str] = None,
    cumulative_path: Optional[str | Path] = None,
    pulser_frequency_hz: Optional[float] = None,
    real_time_sec: Optional[float] = None,
) -> PipelineResult:
    """Convenience: process a single spectrum through the full pipeline.

    See SpectrumPipeline.process() for details.
    """
    pipeline = SpectrumPipeline(
        calibration_path=calibration_path,
        pulser_frequency_hz=pulser_frequency_hz,
        real_time_sec=real_time_sec,
    )
    return pipeline.process(
        raw_path, cumulative_path=cumulative_path, run_id=run_id, source=source
    )


def process_run_sequence(
    run_dir: str | Path,
    calibration_path: str | Path,
    prefix: str = "A",
    run_id: Optional[str] = None,
    source: Optional[str] = None,
    pulser_frequency_hz: Optional[float] = None,
    real_time_sec: Optional[float] = None,
) -> list[PipelineResult]:
    """Convenience: process a run sequence through the full pipeline.

    See SpectrumPipeline.process_sequence() for details.
    """
    pipeline = SpectrumPipeline(
        calibration_path=calibration_path,
        pulser_frequency_hz=pulser_frequency_hz,
        real_time_sec=real_time_sec,
    )
    return pipeline.process_sequence(
        run_dir, prefix=prefix, run_id=run_id, source=source
    )
