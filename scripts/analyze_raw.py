#!/usr/bin/env python3
"""
Analyze raw Tc-99 beta spectra for rate stability and endpoint drift.

Data is cumulative (no reset between runs), so each file contains all
events accumulated up to that point. To get individual ~1-hour runs,
subtract consecutive files: run_N = cumulative_N - cumulative_(N-1).

Produces plots of:
  - Total count rate vs run number
  - Regional count rates (low-energy vs high-energy)
  - Low-energy fraction vs run number
  - Endpoint position vs run number
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import linregress

from exp_data.raw_reader import (
    calibrate_from_am241,
    get_run_files,
    read_raw_spectrum,
    subtract_spectra,
)

# Energy calibration from Am-241 source: E = offset + gain * channel
# Determined automatically from known Am-241 gamma lines (13.8, 17.8, 26.3 keV)
CALIBRATION_FILE = "data/raw/CAM1.DAT"
_calib_offset, _calib_gain = calibrate_from_am241(CALIBRATION_FILE)
ENERGY_OFFSET = _calib_offset
ENERGY_GAIN = _calib_gain
print(f"Energy calibration: E = {ENERGY_OFFSET:.4f} + {ENERGY_GAIN:.4f} * ch")

# Approximate Tc-99 endpoint for reference
TC99_ENDPOINT = 294.0  # keV


def channel_to_energy(ch: np.ndarray) -> np.ndarray:
    """Convert channel index to energy in keV."""
    return ENERGY_OFFSET + ch * ENERGY_GAIN


def estimate_endpoint_from_counts(
    counts: np.ndarray,
    energies_keV: np.ndarray,
    guess_endpoint: float = TC99_ENDPOINT,
) -> float:
    """Estimate beta endpoint from spectrum using Kurie plot linear extrapolation.

    For allowed beta decay, the Kurie plot K(E) = sqrt(counts / (pWF)) is
    linear near the endpoint and extrapolates to zero at E = Q-value.
    We use a simplified Kurie plot (sqrt(counts)) since the phase space
    factor is slowly varying compared to the sharp drop at the endpoint.

    Fits a line to the last ~15% of the spectrum (in the tail region
    where counts are small but above background) and finds the x-intercept.

    Parameters
    ----------
    counts : np.ndarray
        Count spectrum.
    energies_keV : np.ndarray
        Energy calibration for each channel.
    guess_endpoint : float
        Expected endpoint for clamping.

    Returns
    -------
    float
        Estimated endpoint energy in keV.
    """
    # Exclude overflow channel (last channel is always saturated)
    counts_work = counts[:-1]
    energies_work = energies_keV[:-1]

    # Estimate background from the very tail (last 5%)
    n_tail = max(20, int(len(counts_work) * 0.05))
    background = np.median(counts_work[-n_tail:])

    # Select tail region: last 15% of spectrum, above background
    n_tail_region = max(30, int(len(counts_work) * 0.15))
    E_tail = energies_work[-n_tail_region:]
    c_tail = counts_work[-n_tail_region:]

    # Only keep channels above background + noise
    noise = max(np.sqrt(background), 1.0)
    mask = c_tail > (background + noise)
    if not mask.any():
        return guess_endpoint

    E_sel = E_tail[mask]
    c_sel = c_tail[mask]

    # Kurie plot: K = sqrt(counts - background)
    K = np.sqrt(np.maximum(c_sel - background, 0.01))

    # Linear fit: K = a + b * E, endpoint = -a/b
    try:
        slope, intercept, r_val, _, _ = linregress(E_sel, K)[:5]
        if slope < 0 and abs(r_val) > 0.3:
            endpoint = -intercept / slope
            # Clamp to reasonable range (Q ± 100 keV)
            endpoint = max(guess_endpoint - 100, min(guess_endpoint + 100, endpoint))
            return float(endpoint)
    except (ValueError, ZeroDivisionError):
        pass

    # Fallback: use the last channel above background
    last_above = E_sel[np.where(c_sel > background)[0][-1]]
    return float(last_above)


def compute_regional_rates(counts: np.ndarray, energies_keV: np.ndarray) -> dict:
    """Compute count rates in different energy regions."""
    total = float(np.sum(counts))
    low_mask = energies_keV < 100
    high_mask = ~low_mask

    return {
        "total": total,
        "low_energy": float(np.sum(counts[low_mask])),
        "high_energy": float(np.sum(counts[high_mask])),
        "low_fraction": float(np.sum(counts[low_mask]) / total) if total > 0 else 0.0,
    }


def main():
    """Main analysis pipeline."""
    data_dir = Path(__file__).parent.parent / "data" / "raw"
    output_dir = Path(__file__).parent.parent / "output"
    output_dir.mkdir(exist_ok=True)

    # Get run files sorted naturally (A1, A2, ..., A10, A100, ...)
    run_files = get_run_files(data_dir, prefix="A")
    if len(run_files) < 2:
        print(f"Need at least 2 run files, found {len(run_files)}")
        sys.exit(1)

    print(f"Found {len(run_files)} cumulative run files")

    # Step 1: Read all cumulative spectra
    print("\nReading cumulative spectra...")
    cumulative_specs = []
    for i, fpath in enumerate(run_files):
        try:
            spec = read_raw_spectrum(
                fpath,
                run_id=fpath.name,
                source="Tc99",
            )
            cumulative_specs.append(spec)
        except (FileNotFoundError, ValueError) as e:
            print(f"  Skipping {fpath.name}: {e}")
            continue

    if len(cumulative_specs) < 2:
        print("Need at least 2 valid spectra")
        sys.exit(1)

    print(f"  Loaded {len(cumulative_specs)} cumulative spectra")

    # Step 2: Subtract consecutive runs to get individual ~1h runs
    print("\nSubtracting consecutive runs...")
    individual_runs = []
    for i in range(1, len(cumulative_specs)):
        parent_a = cumulative_specs[i].run_id
        parent_b = cumulative_specs[i - 1].run_id
        diff = subtract_spectra(
            cumulative_a=cumulative_specs[i],
            cumulative_b=cumulative_specs[i - 1],
            run_id=f"run_{i}",
        )
        individual_runs.append(diff)

    print(f"  Extracted {len(individual_runs)} individual runs")

    # Step 3: Analyze each run
    print("\nAnalyzing individual runs...")
    run_ids = []
    total_rates = []
    low_rates = []
    high_rates = []
    low_fractions = []
    endpoints = []

    energies_keV = channel_to_energy(np.arange(len(individual_runs[0].counts)))

    for i, run in enumerate(individual_runs):
        run_ids.append(i + 1)
        rates = compute_regional_rates(run.counts, energies_keV)
        endpoint = estimate_endpoint_from_counts(run.counts, energies_keV)

        total_rates.append(rates["total"])
        low_rates.append(rates["low_energy"])
        high_rates.append(rates["high_energy"])
        low_fractions.append(rates["low_fraction"])
        endpoints.append(endpoint)

        if (i + 1) % 50 == 0:
            print(f"  Analyzed {i + 1}/{len(individual_runs)} runs")

    run_ids = np.array(run_ids)
    total_rates = np.array(total_rates)
    low_rates = np.array(low_rates)
    high_rates = np.array(high_rates)
    low_fractions = np.array(low_fractions)
    endpoints = np.array(endpoints)

    # Print summary
    print(f"\n{'='*60}")
    print(f"Analysis of {len(run_ids)} individual Tc-99 runs")
    print(f"{'='*60}")
    print(f"Run range: {run_ids.min()} - {run_ids.max()}")
    print(f"Total counts: {total_rates.min():.0f} - {total_rates.max():.0f}")
    print(f"Low-energy fraction: {low_fractions.min():.3f} - {low_fractions.max():.3f}")
    print(f"Endpoint estimate: {endpoints.min():.1f} - {endpoints.max():.1f} keV")

    # Endpoint drift
    if len(run_ids) > 10:
        slope, intercept, r_val, _, _ = linregress(run_ids, endpoints)
        print(f"\nEndpoint drift: {slope:.3f} keV/run (R²={r_val**2:.4f})")
        print(f"  Start (run 1): {intercept + slope * run_ids[0]:.1f} keV")
        print(f"  End (run {run_ids[-1]}): {intercept + slope * run_ids[-1]:.1f} keV")

    # Rate drift
    slope, intercept, r_val, _, _ = linregress(run_ids, total_rates)
    print(f"\nTotal count rate drift: {slope:.0f} counts/run (R²={r_val**2:.4f})")

    # Save results
    results_file = output_dir / "rate_analysis.csv"
    with open(results_file, "w") as f:
        f.write("run,total_counts,low_energy_counts,high_energy_counts,low_fraction,endpoint_keV\n")
        for i in range(len(run_ids)):
            f.write(
                f"{run_ids[i]},{total_rates[i]:.0f},{low_rates[i]:.0f},"
                f"{high_rates[i]:.0f},{low_fractions[i]:.4f},{endpoints[i]:.1f}\n"
            )
    print(f"\nResults saved to {results_file}")

    # Plot
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Total counts
        axes[0, 0].plot(run_ids, total_rates, "k-", lw=1.5)
        axes[0, 0].set_xlabel("Run number")
        axes[0, 0].set_ylabel("Total counts")
        axes[0, 0].set_title("Total count rate per run")
        axes[0, 0].grid(True, alpha=0.3)

        # Low-energy fraction
        axes[0, 1].plot(run_ids, low_fractions, "b-", lw=1.5)
        axes[0, 1].set_xlabel("Run number")
        axes[0, 1].set_ylabel("Low-energy fraction (<100 keV)")
        axes[0, 1].set_title("Low-energy fraction drift")
        axes[0, 1].grid(True, alpha=0.3)

        # Endpoint
        axes[1, 0].plot(run_ids, endpoints, "r-", lw=1.5)
        axes[1, 0].axhline(
            TC99_ENDPOINT, color="gray", ls="--",
            label=f"Tc-99 Q ≈ {TC99_ENDPOINT} keV",
        )
        axes[1, 0].set_xlabel("Run number")
        axes[1, 0].set_ylabel("Endpoint (keV)")
        axes[1, 0].set_title("Endpoint position drift")
        axes[1, 0].legend(fontsize=9)
        axes[1, 0].grid(True, alpha=0.3)

        # Low vs high energy
        axes[1, 1].plot(run_ids, low_rates, "b-", lw=1.5, label="<100 keV")
        axes[1, 1].plot(run_ids, high_rates, "r-", lw=1.5, label=">100 keV")
        axes[1, 1].set_xlabel("Run number")
        axes[1, 1].set_ylabel("Counts")
        axes[1, 1].set_title("Regional count rates")
        axes[1, 1].legend(fontsize=9)
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plot_file = output_dir / "rate_analysis.png"
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {plot_file}")
    except ImportError:
        print("\nmatplotlib not available — skipping plot generation")


if __name__ == "__main__":
    main()
