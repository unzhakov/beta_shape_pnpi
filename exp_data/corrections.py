"""
Instrumental corrections for experimental spectra.

Provides dead-time correction, pile-up correction, and background
subtraction for beta spectrometry data.
"""

from __future__ import annotations

from typing import Optional

import numpy as np


class DeadTimeCorrection:
    """Dead-time correction for paralyzable or non-paralyzable detectors.

    Parameters
    ----------
    dead_time : float
        Dead time per event in seconds.
    model : str
        'non_paralyzable' (default) or 'paralyzable'.
    """

    def __init__(self, dead_time: float, model: str = "non_paralyzable"):
        self.dead_time = dead_time
        self.model = model
        assert model in (
            "non_paralyzable",
            "paralyzable",
        ), f"Unknown dead-time model: {model}"
        assert dead_time >= 0, "Dead time must be non-negative"

    def live_time_fraction(self, observed_rate: float) -> float:
        """Calculate live time fraction given observed count rate.

        Parameters
        ----------
        observed_rate : float
            Observed count rate in counts/second.

        Returns
        -------
        float
            Fraction of time the detector is live (0 to 1).
        """
        if self.model == "non_paralyzable":
            return max(0.0, 1.0 - observed_rate * self.dead_time)
        else:  # paralyzable
            return np.exp(-observed_rate * self.dead_time)

    def correct_rate(self, observed_rate: float) -> float:
        """Correct observed rate for dead time.

        Parameters
        ----------
        observed_rate : float
            Observed count rate (counts/s).

        Returns
        -------
        float
            True count rate (counts/s).
        """
        live_frac = self.live_time_fraction(observed_rate)
        if live_frac <= 0:
            return float("inf")
        return observed_rate / live_frac

    def correct_counts(self, counts: np.ndarray, live_time: float) -> np.ndarray:
        """Correct counts for dead time.

        Parameters
        ----------
        counts : np.ndarray
            Measured counts per bin.
        live_time : float
            Total live time in seconds.

        Returns
        -------
        np.ndarray
            Dead-time corrected counts.
        """
        total_observed = np.sum(counts)
        total_true = self.correct_rate(total_observed / live_time) * live_time
        if total_observed <= 0:
            return counts.copy()
        return counts * (total_true / total_observed)


class PileUpCorrection:
    """Simple pile-up correction for beta spectrometry.

    Estimates the fraction of pile-up events based on observed count rate
    and dead time, and redistributes pile-up events back into the spectrum.

    This is a first-order approximation. For precise pile-up correction,
    digital pulse processing or anti-pile-up electronics are recommended.
    """

    def __init__(self, dead_time: float):
        self.dead_time = dead_time

    def estimate_pile_up_fraction(self, count_rate: float) -> float:
        """Estimate pile-up fraction from count rate.

        For a non-paralyzable detector, pile-up fraction ≈ (count_rate × dead_time)².

        Parameters
        ----------
        count_rate : float
            Observed count rate in counts/second.

        Returns
        -------
        float
            Estimated fraction of events affected by pile-up.
        """
        return (count_rate * self.dead_time) ** 2

    def correct_spectrum(
        self,
        counts: np.ndarray,
        live_time: float,
        background: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Apply first-order pile-up correction.

        Parameters
        ----------
        counts : np.ndarray
            Measured counts per bin.
        live_time : float
            Live time in seconds.
        background : np.ndarray, optional
            Background counts per bin.

        Returns
        -------
        np.ndarray
            Pile-up corrected counts.
        """
        if background is not None:
            net_counts = counts - background
        else:
            net_counts = counts.copy()

        observed_rate = np.sum(net_counts) / live_time
        pileup_frac = self.estimate_pile_up_fraction(observed_rate)

        # Redistribute pile-up events proportionally to the spectrum shape
        if pileup_frac > 0 and np.sum(net_counts) > 0:
            correction = net_counts * pileup_frac / (1.0 - pileup_frac)
            return net_counts + correction

        return net_counts


class BackgroundSubtractor:
    """Background subtraction for beta spectra.

    Supports:
    - Constant background (from sidebands)
    - Polynomial background
    - Region-of-interest based background estimation
    """

    @staticmethod
    def constant_background(
        counts: np.ndarray,
        energies: np.ndarray,
        low_roi: tuple[float, float],
        high_roi: tuple[float, float],
    ) -> np.ndarray:
        """Subtract constant background estimated from sidebands.

        Parameters
        ----------
        counts : np.ndarray
            Measured counts per bin.
        energies : np.ndarray
            Energy values (keV) per bin.
        low_roi : (E_min, E_max)
            Energy range for low-side background estimation.
        high_roi : (E_min, E_max)
            Energy range for high-side background estimation.

        Returns
        -------
        np.ndarray
            Counts with constant background subtracted.
        """
        low_mask = (energies >= low_roi[0]) & (energies <= low_roi[1])
        high_mask = (energies >= high_roi[0]) & (energies <= high_roi[1])

        low_bg = np.mean(counts[low_mask]) if np.any(low_mask) else 0.0
        high_bg = np.mean(counts[high_mask]) if np.any(high_mask) else 0.0
        constant_bg = (low_bg + high_bg) / 2.0

        return counts - constant_bg

    @staticmethod
    def polynomial_background(
        counts: np.ndarray,
        energies: np.ndarray,
        order: int = 2,
        low_roi: tuple[float, float] = (0.0, 0.0),
        high_roi: tuple[float, float] = (0.0, 0.0),
    ) -> np.ndarray:
        """Subtract polynomial background fitted from sidebands.

        Parameters
        ----------
        counts : np.ndarray
            Measured counts per bin.
        energies : np.ndarray
            Energy values (keV) per bin.
        order : int
            Polynomial order for background fit.
        low_roi : (E_min, E_max)
            Energy range for low-side background fitting.
        high_roi : (E_min, E_max)
            Energy range for high-side background fitting.

        Returns
        -------
        np.ndarray
            Counts with polynomial background subtracted.
        """
        low_mask = (energies >= low_roi[0]) & (energies <= low_roi[1])
        high_mask = (energies >= high_roi[0]) & (energies <= high_roi[1])

        mask = low_mask | high_mask
        if not np.any(mask):
            return counts

        x_fit = energies[mask]
        y_fit = counts[mask]

        coeffs = np.polyfit(x_fit, y_fit, order)
        bg = np.polyval(coeffs[::-1], energies)

        return counts - bg

    @staticmethod
    def roi_average(
        counts: np.ndarray,
        roi_mask: np.ndarray,
    ) -> float:
        """Calculate average background from a region of interest.

        Parameters
        ----------
        counts : np.ndarray
            Measured counts per bin.
        roi_mask : np.ndarray
            Boolean mask selecting background region.

        Returns
        -------
        float
            Average background count per bin.
        """
        if not np.any(roi_mask):
            return 0.0
        return float(np.mean(counts[roi_mask]))
