# spectrum.py
from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from matplotlib.axes import Axes
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple
from dataclasses import dataclass

from beta_spectrum.base import SpectrumComponent
from beta_spectrum.components.phase_space import PhaseSpace
from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.finite_size import FiniteSizeL0, ChargeDistributionU
from beta_spectrum.components.screening import ScreeningCorrection
from beta_spectrum.components.exchange import ExchangeCorrection
from beta_spectrum.components.radiative import RadiativeCorrection

from beta_spectrum.logging_utils import get_git_short_hash
from beta_spectrum.utils import T_to_W
from beta_spectrum.constants import ME_MEV

if TYPE_CHECKING:
    from beta_spectrum.components.detector_response import DetectorResponse


@dataclass
class BranchConfig:
    """
    Configuration for a single decay branch.

    Each branch corresponds to a decay to a specific daughter nuclear state,
    with its own endpoint energy (Q_value - excitation_energy) and transition type.

    Attributes
    ----------
    endpoint_MeV : float
        Effective endpoint energy for this branch (Q_value - excitation_energy).
    transition_type : str
        Forbiddenness classification for this branch (A, F1, F1U, F2, ...).
    intensity : float
        Branch intensity as a fraction of total decays (0.0–1.0).
        Intensities across all branches should sum to 1.0 (100%).
    """

    endpoint_MeV: float
    transition_type: str = "A"
    intensity: float = 1.0

    def __post_init__(self) -> None:
        if self.intensity < 0 or self.intensity > 1.0:
            raise ValueError(
                f"Branch intensity must be in [0, 1], got {self.intensity}"
            )
        if self.endpoint_MeV <= 0:
            raise ValueError(
                f"Branch endpoint_MeV must be positive, got {self.endpoint_MeV}"
            )


@dataclass
class SpectrumConfig:
    """
    Configuration for beta spectrum calculation.

    Supports declarative detector response specification:
    set detector_model, detector_sigma_a_keV, and related parameters
    to enable automatic detector smearing via convolve_with_detector().

    Supports multi-branch decay via the `branches` field. When branches are
    provided, the total spectrum is the intensity-weighted sum of individual
    branch spectra. Each branch has its own endpoint energy and transition type.
    """

    Z_parent: int
    Z_daughter: int
    A_number: int
    endpoint_MeV: float
    transition_type: str = "A"  # should be [A, F1, F1U, F2, F2U, F3, F3U, F4]
    e_step_MeV: float = 0.001

    # Toggle components and corrections
    use_phase_space: bool = True
    use_fermi: bool = True
    use_screening: bool = True
    use_finite_size: bool = True
    use_charge_dist: bool = True
    use_radiative: bool = True
    use_exchange: bool = True

    # Detector response convolution (analytical model)
    use_detector_response: bool = False
    detector_model: str = "gaussian"
    detector_sigma_a_keV: float = 1.0
    detector_sigma_b: float = 0.0
    detector_tail_fraction: float = 0.0
    detector_tau_keV: float = 5.0
    detector_fano_factor: float = 0.12
    detector_n_channels: int = 4096
    detector_channel_energy_range: tuple[float, float] = (
        0.0,
        0.35,
    )  # in m_e units (total energy)

    # Multi-branch support
    branches: Optional[List[BranchConfig]] = None  # None → single-branch mode
    intensity_cutoff: float = (
        0.0  # Ignore branches below this fraction (default: 0.0 = all)
    )


class BetaSpectrum:
    """
    Precise beta spectrum calculator.

    Supports both single-branch and multi-branch decay modes.

    In single-branch mode (branches=None), the spectrum is the product of
    all enabled correction factors.

    In multi-branch mode (branches provided), each branch is calculated
    independently with its own endpoint energy and transition type, then
    summed with intensity weighting.

    Parameters
    ----------
    components : List[SpectrumComponent], optional
        List of enabled universal spectral components (Fermi, screening, etc.).
        Only used in multi-branch mode or for universal components.
    branch_spectra : List[BetaSpectrum], optional
        List of per-branch spectrum calculators (each with its own PhaseSpace, etc.).
        Only used in multi-branch mode.
    branches : List[BranchConfig], optional
        Branch configuration metadata (intensities, endpoint energies).
    logger : logging.Logger, optional
        Logger for debug/info output.
    """

    def __init__(
        self,
        components: Optional[List[SpectrumComponent]] = None,
        branch_spectra: Optional[List[BetaSpectrum]] = None,
        branches: Optional[List[BranchConfig]] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.components = components or []
        self._component_names = [self._get_component_name(c) for c in self.components]
        self.branch_spectra = branch_spectra or []
        self.branches = branches or []
        self._logger = logger
        self._is_multi = len(self.branch_spectra) > 0

    def __call__(self, W: np.ndarray) -> np.ndarray:
        """
        Calculate spectrum.

        In single-branch mode: product of all components.
        In multi-branch mode: intensity-weighted sum of branch spectra.
        """
        if self._is_multi:
            return self._multi_branch_call(W)
        return self._single_branch_call(W)

    def _single_branch_call(self, W: np.ndarray) -> np.ndarray:
        """Calculate spectrum as product of components (single branch)."""
        result: np.ndarray = np.ones_like(W, dtype=float)
        for comp in self.components:
            if self._logger:
                self._logger.debug(
                    "Evaluating %s at %d energy points", comp.__class__.__name__, len(W)
                )
            result = result * comp(W)
        return result.astype(np.float64)

    def _multi_branch_call(self, W: np.ndarray) -> np.ndarray:
        """Calculate total spectrum as intensity-weighted sum of branches.

        The total is: (universal_components_product) * sum(intensity_i * branch_spectrum_i)
        """
        # Sum of intensity-weighted branch spectra
        branch_sum = np.zeros_like(W, dtype=float)
        for i, branch_spec in enumerate(self.branch_spectra):
            if self._logger:
                self._logger.debug(
                    "Evaluating branch %d/%d (endpoint=%.3f MeV, intensity=%.4f)",
                    i + 1,
                    len(self.branch_spectra),
                    self.branches[i].endpoint_MeV,
                    self.branches[i].intensity,
                )
            branch_spectrum = branch_spec(W)
            # Mask out energies beyond this branch's endpoint
            W0_branch = T_to_W(self.branches[i].endpoint_MeV)
            mask = W <= W0_branch
            branch_sum += self.branches[i].intensity * np.where(
                mask, branch_spectrum, 0.0
            )

        # Multiply by universal components
        result: np.ndarray = np.ones_like(W, dtype=float)
        for comp in self.components:
            result = result * comp(W)

        return np.asarray((result * branch_sum).astype(np.float64), dtype=np.float64)

    def _get_component_name(self, comp: SpectrumComponent) -> str:
        """Get human-readable component name."""
        class_name = comp.__class__.__name__
        return class_name.replace("Correction", "").replace("Function", "")

    def calculate_components(self, W: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate individual components.

        In single-branch mode: returns per-component values.
        In multi-branch mode: returns universal components + per-branch
        branch-dependent components (PhaseSpace, Radiative).
        """
        if self._is_multi:
            return self._multi_branch_components(W)

        result: Dict[str, np.ndarray] = {}
        for name, comp in zip(self._component_names, self.components):
            result[name] = comp(W)
        return result

    def _multi_branch_components(self, W: np.ndarray) -> Dict[str, np.ndarray]:
        """Return universal + per-branch components for multi-branch mode."""
        result: Dict[str, np.ndarray] = {}

        # Universal components (same for all branches)
        for name, comp in zip(self._component_names, self.components):
            result[name] = comp(W)

        # Per-branch components (branch-dependent)
        for i, branch_spec in enumerate(self.branch_spectra):
            branch_name = f"branch_{i}"
            branch_components = branch_spec.calculate_components(W)
            for comp_name, values in branch_components.items():
                result[f"{branch_name}.{comp_name}"] = values

        return result

    @classmethod
    def from_config(
        cls, config: SpectrumConfig, logger: Optional[logging.Logger] = None
    ) -> BetaSpectrum:
        """
        Create a BetaSpectrum from configuration.

        In single-branch mode: creates a single spectrum calculator.
        In multi-branch mode: creates per-branch calculators + universal components.

        Parameters
        ----------
        config : SpectrumConfig
            Configuration with decay parameters and correction toggles.
        logger : logging.Logger, optional
            Logger for progress output.
        """
        if config.branches is not None and len(config.branches) > 1:
            return cls._from_config_multi(config, logger)
        return cls._from_config_single(config, logger)

    @classmethod
    def _from_config_single(
        cls, config: SpectrumConfig, logger: Optional[logging.Logger] = None
    ) -> BetaSpectrum:
        """Create single-branch spectrum calculator."""
        components: List[SpectrumComponent] = []
        W0 = float(T_to_W(config.endpoint_MeV))

        if config.use_phase_space:
            components.append(
                PhaseSpace(W0=W0, transition_type=config.transition_type, logger=logger)
            )

        if config.use_fermi:
            components.append(
                FermiFunction(Z=config.Z_daughter, A=config.A_number, logger=logger)
            )

        if config.use_finite_size:
            components.append(
                FiniteSizeL0(Z=config.Z_daughter, A=config.A_number, logger=logger)
            )

        if config.use_charge_dist:
            components.append(
                ChargeDistributionU(
                    Z=config.Z_daughter, A=config.A_number, logger=logger
                )
            )

        if config.use_screening:
            components.append(
                ScreeningCorrection(
                    FermiFunction(Z=config.Z_parent, A=config.A_number, logger=logger),
                    logger=logger,
                )
            )

        if config.use_exchange:
            components.append(ExchangeCorrection(Z=config.Z_parent, logger=logger))

        if config.use_radiative:
            components.append(
                RadiativeCorrection(
                    W0=W0,
                    Z=config.Z_daughter,
                    A=config.A_number,
                    use_endpoint_resummation=True,
                    logger=logger,
                )
            )

        return cls(components, logger=logger)

    @classmethod
    def _from_config_multi(
        cls, config: SpectrumConfig, logger: Optional[logging.Logger] = None
    ) -> BetaSpectrum:
        """Create multi-branch spectrum calculator.

        Creates per-branch calculators (each with its own PhaseSpace, Radiative, etc.)
        and universal components (Fermi, screening, finite_size, etc.) shared across
        all branches.

        Note: Radiative correction is branch-dependent (depends on W0), so it goes
        into each branch's calculator, not the universal components.
        """
        if logger:
            logger.info("Multi-branch mode: %d branches", len(config.branches))

        # Create per-branch calculators
        branch_spectra: List[BetaSpectrum] = []
        for i, branch in enumerate(config.branches):
            if logger:
                logger.info(
                    "  Branch %d: endpoint=%.3f MeV, transition=%s, intensity=%.4f",
                    i + 1,
                    branch.endpoint_MeV,
                    branch.transition_type,
                    branch.intensity,
                )

            # Create a per-branch config with the branch's endpoint and transition type
            branch_config = SpectrumConfig(
                Z_parent=config.Z_parent,
                Z_daughter=config.Z_daughter,
                A_number=config.A_number,
                endpoint_MeV=branch.endpoint_MeV,
                transition_type=branch.transition_type,
                e_step_MeV=config.e_step_MeV,
                # Branch-dependent components
                use_phase_space=config.use_phase_space,
                use_radiative=config.use_radiative,
                # Universal components disabled (handled separately)
                use_fermi=False,
                use_screening=False,
                use_finite_size=False,
                use_charge_dist=False,
                use_exchange=False,
            )

            branch_spectrum = cls._from_config_single(branch_config, logger)
            branch_spectra.append(branch_spectrum)

        # Create universal components (shared across all branches)
        # Note: Radiative is NOT included here because it depends on W0
        universal_components: List[SpectrumComponent] = []

        if config.use_fermi:
            universal_components.append(
                FermiFunction(Z=config.Z_daughter, A=config.A_number, logger=logger)
            )

        if config.use_finite_size:
            universal_components.append(
                FiniteSizeL0(Z=config.Z_daughter, A=config.A_number, logger=logger)
            )

        if config.use_charge_dist:
            universal_components.append(
                ChargeDistributionU(
                    Z=config.Z_daughter, A=config.A_number, logger=logger
                )
            )

        if config.use_screening:
            universal_components.append(
                ScreeningCorrection(
                    FermiFunction(Z=config.Z_parent, A=config.A_number, logger=logger),
                    logger=logger,
                )
            )

        if config.use_exchange:
            universal_components.append(
                ExchangeCorrection(Z=config.Z_parent, logger=logger)
            )

        return cls(
            components=universal_components,
            branch_spectra=branch_spectra,
            branches=list(config.branches),
            logger=logger,
        )

    def get_energy_grid(self, config: SpectrumConfig) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate array of energy values for spectrum calculation.

        In multi-branch mode, uses the maximum endpoint across all branches.
        """
        if self._is_multi and self.branches:
            max_endpoint = max(b.endpoint_MeV for b in self.branches)
        else:
            max_endpoint = config.endpoint_MeV

        kinetic_MeV = np.arange(config.e_step_MeV, max_endpoint, config.e_step_MeV)
        W = T_to_W(kinetic_MeV)

        return W, kinetic_MeV.astype(np.float64)

    def get_branch_spectra(self, W: np.ndarray) -> List[np.ndarray]:
        """
        Get individual branch spectra for multi-branch mode.

        Each branch spectrum is masked to zero beyond its endpoint.

        Parameters
        ----------
        W : np.ndarray
            Energy grid in m_e units.

        Returns
        -------
        List[np.ndarray]
            Per-branch spectrum arrays (unweighted). Empty in single-branch mode.
        """
        if not self._is_multi:
            return []
        result = []
        for i, branch_spec in enumerate(self.branch_spectra):
            branch_spectrum = branch_spec(W)
            # Use np.maximum to smoothly approach zero at endpoint,
            # and clamp negative values beyond endpoint to zero
            result.append(np.maximum(0.0, branch_spectrum))
        return result

    def get_branch_normalized_spectra(self, W: np.ndarray) -> List[np.ndarray]:
        """
        Get intensity-weighted branch spectra for multi-branch mode.

        Each branch spectrum is masked to zero beyond its endpoint.

        Parameters
        ----------
        W : np.ndarray
            Energy grid in m_e units.

        Returns
        -------
        List[np.ndarray]
            Per-branch spectrum arrays weighted by intensity. Empty in single-branch mode.
        """
        if not self._is_multi:
            return []
        result = []
        for i, branch_spec in enumerate(self.branch_spectra):
            branch_spectrum = branch_spec(W)
            result.append(self.branches[i].intensity * np.maximum(0.0, branch_spectrum))
        return result

    @staticmethod
    def create_detector_from_config(
        config: SpectrumConfig,
    ) -> "DetectorResponse":
        """
        Create a DetectorResponse from SpectrumConfig detector parameters.

        Converts keV-based resolution parameters to m_e units internally.

        Parameters
        ----------
        config : SpectrumConfig
            Configuration with detector response parameters set.

        Returns
        -------
        DetectorResponse
            Detector response object ready for convolution.
        """
        from beta_spectrum.components.detector_response import DetectorResponse
        from beta_spectrum.utils import T_to_W

        sigma_a_me = config.detector_sigma_a_keV / ME_MEV

        tau_me = config.detector_tau_keV / ME_MEV

        W0 = T_to_W(config.endpoint_MeV)
        channel_range = config.detector_channel_energy_range

        # Extend range to cover endpoint
        if channel_range[1] < W0:
            channel_range = (channel_range[0], float(W0 + 0.05))

        detector = DetectorResponse.from_gaussian_params(
            channel_energy_range=channel_range,
            n_channels=config.detector_n_channels,
            sigma_a=sigma_a_me,
            sigma_b=config.detector_sigma_b,
            tail_fraction=config.detector_tail_fraction,
            tau=tau_me,
            model=config.detector_model,
            fano_factor=config.detector_fano_factor,
        )
        return detector

    def convolve_with_detector(
        self,
        detector_response: "DetectorResponse",
        W: Optional[np.ndarray] = None,
        config: Optional[SpectrumConfig] = None,
    ) -> np.ndarray:
        """
        Convolve theoretical spectrum with detector response.

        Returns the predicted measured spectrum after detector smearing.

        Parameters
        ----------
        detector_response : DetectorResponse
            Detector response object (analytical or tabulated).
        W : np.ndarray, optional
            Energy grid in m_e units. If None, generated from config.
        config : SpectrumConfig, optional
            Used to generate energy grid if W is None.

        Returns
        -------
        np.ndarray
            Convolved spectrum (predicted measured counts per channel).
        """
        if W is None:
            if config is None:
                raise ValueError("Either W or config must be provided")
            W, _ = self.get_energy_grid(config)

        theoretical_spectrum = self(W)

        if self._logger:
            self._logger.debug(
                "Convolution: %d channels, %d spectrum points",
                len(detector_response.channel_energies),
                len(W),
            )

        convolved = detector_response.convolve(W, theoretical_spectrum, normalize=True)
        return convolved

    def convolve_detector(
        self,
        config: SpectrumConfig,
        W: Optional[np.ndarray] = None,
        detector_response: Optional["DetectorResponse"] = None,
    ) -> np.ndarray:
        """
        Convolve theoretical spectrum with detector response from config.

        Convenience method: creates DetectorResponse from config parameters
        and convolves the spectrum in a single call.

        Parameters
        ----------
        config : SpectrumConfig
            Configuration with detector response parameters set.
        W : np.ndarray, optional
            Energy grid in m_e units. If None, generated from config.
        detector_response : DetectorResponse, optional
            Pre-built detector response. If None, created from config.

        Returns
        -------
        np.ndarray
            Convolved spectrum (predicted measured counts per channel).

        Examples
        --------
        >>> config = SpectrumConfig(
        ...     Z_parent=43, Z_daughter=44, A_number=99, endpoint_MeV=0.294,
        ...     use_detector_response=True,
        ...     detector_sigma_a_keV=1.0,
        ... )
        >>> spectrum = BetaSpectrum.from_config(config)
        >>> W, _ = spectrum.get_energy_grid(config)
        >>> convolved = spectrum.convolve_detector(config, W=W)
        """
        if detector_response is None:
            detector_response = self.create_detector_from_config(config)

        return self.convolve_with_detector(detector_response, W=W, config=config)


class BetaSpectrumAnalyzer:
    """
    Introspection, analysis and debugging tools for BetaSpectrum.

    Parameters
    ----------
    spectrum : BetaSpectrum
        The spectrum to analyze.
    config : SpectrumConfig
        Configuration used to create the spectrum.
    logger : logging.Logger, optional
        Logger for progress output.
    """

    def __init__(
        self,
        spectrum: BetaSpectrum,
        config: SpectrumConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Initialize analyzer with a spectrum and configuration.
        """
        self.spectrum = spectrum
        self.config = config
        self.W, self.energies_MeV = spectrum.get_energy_grid(config)
        self._components_cache: Optional[Dict[str, np.ndarray]] = None
        self._branch_spectra_cache: Optional[List[np.ndarray]] = None
        self._logger = logger
        self._is_multi = spectrum._is_multi

    @property
    def components(self) -> Dict[str, np.ndarray]:
        if self._components_cache is None:
            self._components_cache = self.spectrum.calculate_components(self.W)
        assert self._components_cache is not None
        return self._components_cache

    @property
    def branch_spectra(self) -> List[np.ndarray]:
        """Intensity-weighted branch spectra."""
        if self._branch_spectra_cache is None:
            self._branch_spectra_cache = self.spectrum.get_branch_normalized_spectra(
                self.W
            )
        assert self._branch_spectra_cache is not None
        return self._branch_spectra_cache

    def total_spectrum(self, normalize: bool = True) -> np.ndarray:
        """
        Calculate the total_spectrum
        """
        total = self.spectrum(self.W)

        if normalize:
            integral = np.trapezoid(total, self.energies_MeV)
            total = total / integral

        return total

    def convolved_spectrum(
        self,
        detector_response: Optional["DetectorResponse"] = None,
        normalize: bool = True,
    ) -> np.ndarray:
        """
        Calculate spectrum convolved with detector response.

        If detector_response is None and use_detector_response is True in config,
        creates a detector response from config parameters.

        Parameters
        ----------
        detector_response : DetectorResponse, optional
            Detector response to convolve with. Created from config if None.
        normalize : bool
            If True, normalize to unit area.

        Returns
        -------
        np.ndarray
            Convolved spectrum.
        """
        if detector_response is None and self.config.use_detector_response:
            detector_response = BetaSpectrum.create_detector_from_config(self.config)

        if detector_response is None:
            raise ValueError(
                "detector_response must be provided, or set "
                "use_detector_response=True in config"
            )

        convolved = self.spectrum.convolve_with_detector(
            detector_response, W=self.W, config=self.config
        )

        if normalize:
            integral = np.trapezoid(convolved, detector_response.channel_energies)
            if integral > 0:
                convolved = convolved / integral

        return convolved

    def plot_analysis(
        self,
        save_path: Optional[str] = None,
        show_components: bool = True,
    ) -> None:
        """
        Create visualization of the spectrum and all correction factors.

        In single-branch mode: standard 4-panel debug view or single spectrum plot.
        In multi-branch mode: spectrum plot with individual branches shown,
        and debug view with vertical layout for branch comparison.

        Parameters
        ----------
        save_path : str, optional
            Path to save the figure.
        show_components : bool
            If True, show full debug view.
            If False, show only the spectrum plot with nuclear data header.
            If True, use logarithmic y-scale. Set to False for
            linear scale, which is better when branch intensities are similar.
        """
        total = self.total_spectrum(normalize=True)
        components = self.components
        commit = get_git_short_hash(6)
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        if self._is_multi:
            if show_components:
                self._plot_multi_branch_debug(
                    total, components, save_path, commit, timestamp
                )
            else:
                self._plot_multi_branch_spectra(
                    total, components, save_path, commit, timestamp
                )
        else:
            if show_components:
                self._plot_debug_view(total, components, save_path, commit, timestamp)
            else:
                self._plot_spectrum_only(total, save_path, commit, timestamp)

    def _add_id_textbox(self, ax: "Axes", commit: str, timestamp: str) -> None:
        """Add ID text box with commit hash and timestamp to plot."""
        id_text = f"commit: {commit}  |  {timestamp}"
        ax.text(
            0.01,
            0.01,
            id_text,
            transform=ax.transAxes,
            fontsize=8,
            color="gray",
            verticalalignment="bottom",
            horizontalalignment="left",
            fontfamily="monospace",
        )

    def _add_nuclear_data_header(self, ax: "Axes") -> None:
        """Add nuclear data information header to plot."""
        parent = self._element_symbol(self.config.Z_parent)
        daughter = self._element_symbol(self.config.Z_daughter)

        enabled_corrections = []
        if self.config.use_phase_space:
            enabled_corrections.append("phase_space")
        if self.config.use_fermi:
            enabled_corrections.append("fermi")
        if self.config.use_screening:
            enabled_corrections.append("screening")
        if self.config.use_finite_size:
            enabled_corrections.append("finite_size")
        if self.config.use_charge_dist:
            enabled_corrections.append("charge_dist")
        if self.config.use_radiative:
            enabled_corrections.append("radiative")
        if self.config.use_exchange:
            enabled_corrections.append("exchange")

        header_lines = [
            f"Nuclide:    {parent}{self.config.A_number} -> {daughter}{self.config.A_number}",
            f"Endpoint:   {self.config.endpoint_MeV * 1000:.1f} keV",
            f"Transition: {self.config.transition_type}",
            f"Corrections: {', '.join(enabled_corrections)}",
        ]

        if self.config.use_detector_response:
            header_lines.append(
                f"Detector:   {self.config.detector_model} "
                f"(σ={self.config.detector_sigma_a_keV} keV)"
            )

        header_text = "\n".join(header_lines)
        ax.text(
            0.99,
            0.99,
            header_text,
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"
            ),
        )

    def _plot_spectrum_only(
        self,
        total: np.ndarray,
        save_path: Optional[str],
        commit: str,
        timestamp: str,
    ) -> None:
        """Plot only the total spectrum with nuclear data header."""
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(self.energies_MeV, total, "b-", lw=2, label="Normalized spectrum")
        ax.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=11)
        ax.set_ylabel("Normalized Counts", fontsize=11)
        ax.set_title(
            f"Beta-decay: {self.config.Z_parent} -> {self.config.Z_daughter}, A={self.config.A_number}",
            fontsize=13,
        )
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=10)

        self._add_nuclear_data_header(ax)
        self._add_id_textbox(ax, commit, timestamp)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

        plt.show()

    def _plot_debug_view(
        self,
        total: np.ndarray,
        components: Dict[str, np.ndarray],
        save_path: Optional[str],
        commit: str,
        timestamp: str,
    ) -> None:
        """Full debug view with 4-panel spectrum analysis."""
        fig = plt.figure(figsize=(14, 10))

        # 1. Main spectrum plot (top-left)
        ax1 = plt.subplot(2, 2, 1)
        ax1.plot(
            self.energies_MeV, total, "b-", lw=2, label="Total spectrum (log scale)"
        )
        ax1.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=10)
        ax1.set_ylabel(
            f"Normalized Counts per {self.config.e_step_MeV:.3e} MeV", fontsize=10
        )
        ax1.set_title(
            f"Beta-decay: Z={self.config.Z_parent} -> {self.config.Z_daughter}, A={self.config.A_number}",
            fontsize=12,
        )
        ax1.set_yscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # 2. Correction plots (top-right)
        ax2 = plt.subplot(2, 2, 2)

        # Set color order for correction plots
        component_order = [
            "PhaseSpace",
            "Fermi",
            "Screening",
            "Exchange",
            "FiniteSizeL0",
            "ChargeDistributionU",
            "Radiative",
        ]
        colors = ["gray", "red", "green", "blue", "orange", "purple", "brown"]

        for name, color in zip(component_order, colors):
            if name in components:
                values = components[name]

                # For Fermi and PhaseSpace, scale for visibility
                if name in ["Fermi", "PhaseSpace"]:
                    values = values / np.max(values)
                    label = f"{name} (norm)"
                else:
                    label = name

                ax2.plot(
                    self.energies_MeV,
                    values,
                    color=color,
                    lw=1.5,
                    label=label,
                    alpha=0.8,
                )

        ax2.set_xlabel("Electron kinetic energy E [MeV]", fontsize=10)
        ax2.set_ylabel("Correction factor", fontsize=10)
        ax2.set_title("Spectrum components", fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc="best", fontsize=8)

        # 3. Cumulative effect (botom-left)
        ax3 = plt.subplot(2, 2, 3)

        # Cumulative product
        cumulative = np.ones_like(self.W)
        baseline_norm = np.ones_like(cumulative)

        # Plot baseline
        ax3.plot(
            self.energies_MeV, baseline_norm, "k--", lw=1.5, label="Baseline", alpha=0.7
        )

        # Add components one by one
        for name, color in zip(component_order, colors):
            if name in components:
                cumulative *= components[name]
                norm_cumulative = cumulative / np.max(cumulative)
                ax3.plot(
                    self.energies_MeV,
                    norm_cumulative,
                    color=color,
                    lw=1.5,
                    label=f"+ {name}",
                    alpha=0.7,
                )

        ax3.set_xlabel("Electron kinetic energy E [MeV]", fontsize=10)
        ax3.set_ylabel("Normalized spectrum", fontsize=10)
        ax3.set_title("Cumulative effect", fontsize=12)
        ax3.grid(True, alpha=0.3)
        ax3.legend(loc="best", fontsize=8)

        # 4. Deviation from unity (bottom-right)
        ax4 = plt.subplot(2, 2, 4)

        for name, color in zip(component_order, colors):
            if name in components and name not in ["Fermi", "PhaseSpace"]:
                deviation = components[name] - 1.0
                ax4.plot(
                    self.energies_MeV,
                    deviation,
                    color=color,
                    lw=1.5,
                    label=f"{name}",
                    alpha=0.7,
                )

        ax4.set_xlabel("Electron kinetic energy E [MeV]", fontsize=10)
        ax4.set_ylabel("Deviation from unity", fontsize=10)
        ax4.set_title("Correction deviations (C - 1)", fontsize=12)
        ax4.grid(True, alpha=0.3)
        ax4.legend(loc="best", fontsize=8)
        ax4.axhline(y=0, color="k", linestyle="-", lw=0.5)

        plt.suptitle(
            f"Beta-decay spectrum {self.config.A_number}{self._element_symbol(self.config.Z_parent)} -> {self.config.A_number}{self._element_symbol(self.config.Z_daughter)}",
            fontsize=14,
            fontweight="bold",
        )
        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

        plt.show()

    def _element_symbol(self, Z: int) -> str:
        symbols = [
            "",
            "H",
            "He",
            "Li",
            "Be",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
            "K",
            "Ca",
            "Sc",
            "Ti",
            "V",
            "Cr",
            "Mn",
            "Fe",
            "Co",
            "Ni",
            "Cu",
            "Zn",
            "Ga",
            "Ge",
            "As",
            "Se",
            "Br",
            "Kr",
            "Rb",
            "Sr",
            "Y",
            "Zr",
            "Nb",
            "Mo",
            "Tc",
            "Ru",
            "Rh",
            "Pd",
            "Ag",
            "Cd",
            "In",
            "Sn",
            "Sb",
            "Te",
            "I",
            "Xe",
            "Cs",
            "Ba",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "Hf",
            "Ta",
            "W",
            "Re",
            "Os",
            "Ir",
            "Pt",
            "Au",
            "Hg",
            "Tl",
            "Pb",
            "Bi",
            "Po",
            "At",
            "Rn",
            "Fr",
            "Ra",
            "Ac",
            "Th",
            "Pa",
            "U",
            "Np",
            "Pu",
            "Am",
            "Cm",
            "Bk",
            "Cf",
            "Es",
            "Fm",
        ]
        return symbols[Z] if Z < len(symbols) else f"Z{Z}"

    # ---------------------------------------------------------------------------
    # Multi-branch plotting
    # ---------------------------------------------------------------------------

    def _get_branch_label(self, i: int) -> str:
        """Generate branch label for legend."""
        if i >= len(self.spectrum.branches):
            return f"Branch {i+1}"
        branch = self.spectrum.branches[i]
        return f"Branch {i+1}: {branch.transition_type}, E₀={branch.endpoint_MeV*1000:.1f} keV"

    def _plot_multi_branch_spectra(
        self,
        total: np.ndarray,
        components: Dict[str, np.ndarray],
        save_path: Optional[str],
        commit: str,
        timestamp: str,
    ) -> None:
        """Plot total spectrum with individual branch decomposition.

        All curves use raw (un-normalized) values so that the total is
        the arithmetic sum of branch contributions, and branch amplitudes
        are proportional to their intensities.
        """
        fig, ax = plt.subplots(figsize=(12, 7))

        branch_colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]

        # Compute branch contributions (universal * intensity * branch_spectrum)
        branch_contribs = []
        for i, branch_spec in enumerate(self.spectrum.branch_spectra):
            universal_product = np.ones_like(self.W)
            for name in [
                "Fermi",
                "Screening",
                "FiniteSizeL0",
                "ChargeDistributionU",
                "Exchange",
            ]:
                if name in components:
                    universal_product *= components[name]
            bs = self.spectrum.get_branch_spectra(self.W)
            contrib = universal_product * self.spectrum.branches[i].intensity * bs[i]
            branch_contribs.append(contrib)

        # Total = arithmetic sum of all branches (raw, un-normalized)
        total_raw = sum(branch_contribs)

        ax.plot(
            self.energies_MeV,
            total_raw,
            "k-",
            lw=2.5,
            label="Total spectrum",
            zorder=10,
        )

        # Plot each branch contribution (raw, un-normalized)
        for i, (contrib, branch) in enumerate(
            zip(branch_contribs, self.spectrum.branches)
        ):
            color = branch_colors[i % len(branch_colors)]
            ax.plot(
                self.energies_MeV,
                contrib,
                color=color,
                lw=1.5,
                alpha=0.7,
                label=self._get_branch_label(i),
            )

        ax.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=12)
        ax.set_ylabel("Normalized Counts", fontsize=12)
        parent = self._element_symbol(self.config.Z_parent)
        daughter = self._element_symbol(self.config.Z_daughter)
        ax.set_title(
            f"Beta-decay: {parent}{self.config.A_number} -> {daughter}{self.config.A_number} — {len(self.spectrum.branches)} branches",
            fontsize=14,
        )
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=10, ncol=min(2, len(self.spectrum.branches) + 1))

        self._add_nuclear_data_header(ax)
        self._add_id_textbox(ax, commit, timestamp)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

    def _plot_multi_branch_debug(
        self,
        total: np.ndarray,
        components: Dict[str, np.ndarray],
        save_path: Optional[str],
        commit: str,
        timestamp: str,
    ) -> None:
        """Debug view for multi-branch: vertical layout with 4 panels.

        Panel 1: Total spectrum + intensity-scaled branch spectra (log scale)
        Panel 2: Fermi function + all phase spaces (shape comparison, no intensity)
        Panel 3: Universal W0-independent corrections (Fermi, Screening, etc.)
        Panel 4: W0-dependent corrections (Radiative for each branch)
        """
        n_branches = len(self.spectrum.branches)

        fig, axes = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

        branch_colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        parent = self._element_symbol(self.config.Z_parent)
        daughter = self._element_symbol(self.config.Z_daughter)

        # Panel 1: Total spectrum + per-branch spectra
        # All curves use raw (un-normalized) values so total = sum of branches
        branch_contribs = []
        for i, branch_spec in enumerate(self.spectrum.branch_spectra):
            universal_product = np.ones_like(self.W)
            for name in [
                "Fermi",
                "Screening",
                "FiniteSizeL0",
                "ChargeDistributionU",
                "Exchange",
            ]:
                if name in components:
                    universal_product *= components[name]
            bs = self.spectrum.get_branch_spectra(self.W)
            contrib = universal_product * self.spectrum.branches[i].intensity * bs[i]
            branch_contribs.append(contrib)

        total_raw = sum(branch_contribs)

        ax = axes[0]
        ax.plot(
            self.energies_MeV,
            total_raw,
            "k-",
            lw=2,
            label="Total spectrum",
            zorder=10,
        )

        for i, (contrib, branch) in enumerate(
            zip(branch_contribs, self.spectrum.branches)
        ):
            color = branch_colors[i % len(branch_colors)]
            ax.plot(
                self.energies_MeV,
                contrib,
                color=color,
                lw=1.5,
                alpha=0.7,
                label=f"Branch {i+1}: {branch.transition_type}, "
                f"E₀={branch.endpoint_MeV*1000:.1f} keV",
            )
        ax.set_yscale("log")
        ax.set_ylabel("Normalized Counts", fontsize=10)
        ax.set_title(
            f"Multi-branch: {parent}{self.config.A_number} -> {daughter}{self.config.A_number} ({n_branches} branches)",
            fontsize=12,
            fontweight="bold",
        )
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best", ncol=min(2, n_branches + 1))

        # Panel 2: Fermi function + all phase spaces (shape comparison)
        # Phase space is masked beyond endpoint before normalizing to avoid
        # negative values from the sqrt term
        ax = axes[1]
        if "Fermi" in components:
            fermi_vals = components["Fermi"]
            fermi_scaled = fermi_vals / np.max(fermi_vals)
            ax.plot(
                self.energies_MeV,
                fermi_scaled,
                "r-",
                lw=1.5,
                label="Fermi (normalized)",
                alpha=0.8,
            )
        for i, (branch_spec, branch) in enumerate(
            zip(self.spectrum.branch_spectra, self.spectrum.branches)
        ):
            color = branch_colors[i % len(branch_colors)]
            if f"branch_{i}.PhaseSpace" in components:
                ps_vals = components[f"branch_{i}.PhaseSpace"]
                # Mask beyond endpoint to avoid negative values
                W0_branch = T_to_W(branch.endpoint_MeV)
                mask = self.W <= W0_branch
                ps_masked = np.where(mask, ps_vals, 0.0)
                ps_max = np.max(ps_masked)
                if ps_max > 0:
                    ps_scaled = ps_masked / ps_max
                else:
                    ps_scaled = ps_masked
                ax.plot(
                    self.energies_MeV,
                    ps_scaled,
                    color=color,
                    lw=1.5,
                    alpha=0.7,
                    label=f"PS Branch {i+1} (normalized)",
                )
        ax.set_ylabel("Normalized Factor", fontsize=10)
        ax.set_title(
            "Fermi Function & Phase Space Shapes", fontsize=11, fontweight="bold"
        )
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")

        # Panel 3: Universal W0-independent corrections
        ax = axes[2]
        universal_components = [
            "Screening",
            "FiniteSizeL0",
            "ChargeDistributionU",
            "Exchange",
        ]
        colors = ["green", "orange", "purple", "brown"]
        for name, color in zip(universal_components, colors):
            if name in components:
                values = components[name]
                values_scaled = values / np.max(values)
                ax.plot(
                    self.energies_MeV,
                    values_scaled,
                    color=color,
                    lw=1.5,
                    label=f"{name} (normalized)",
                    alpha=0.8,
                )
        ax.set_ylabel("Normalized Factor", fontsize=10)
        ax.set_title(
            "Universal W0-Independent Corrections", fontsize=11, fontweight="bold"
        )
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")

        # Panel 4: W0-dependent corrections (Radiative for each branch)
        ax = axes[3]
        for i, (branch_spec, branch) in enumerate(
            zip(self.spectrum.branch_spectra, self.spectrum.branches)
        ):
            color = branch_colors[i % len(branch_colors)]
            if f"branch_{i}.Radiative" in components:
                rad_vals = components[f"branch_{i}.Radiative"]
                ax.plot(
                    self.energies_MeV,
                    rad_vals,
                    color=color,
                    lw=1.5,
                    alpha=0.7,
                    label=f"Branch {i+1} Radiative (E₀={branch.endpoint_MeV*1000:.1f} keV)",
                )
        ax.set_ylabel("Radiative Correction", fontsize=10)
        ax.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=10)
        ax.set_title(
            "W0-Dependent Corrections (Radiative)", fontsize=11, fontweight="bold"
        )
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best", ncol=min(2, n_branches + 1))

        self._add_id_textbox(axes[3], commit, timestamp)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

    def export_to_csv(self, filename: str, source_type: str = "unknown") -> None:
        """
        Export spectrum data to CSV with metadata header.

        Parameters
        ----------
        filename : str
            Output CSV file path.
        source_type : str
            Data source type (paceENSDF, json, cli) for the header.
        """
        total = self.total_spectrum(normalize=True)
        components = self.components

        # Build enabled corrections list
        enabled = []
        if self.config.use_phase_space:
            enabled.append("phase_space")
        if self.config.use_fermi:
            enabled.append("fermi")
        if self.config.use_screening:
            enabled.append("screening")
        if self.config.use_finite_size:
            enabled.append("finite_size")
        if self.config.use_charge_dist:
            enabled.append("charge_dist")
        if self.config.use_radiative:
            enabled.append("radiative")
        if self.config.use_exchange:
            enabled.append("exchange")

        # Build metadata header
        from datetime import datetime, timezone

        def _element_symbol(Z: int) -> str:
            symbols = [
                "",
                "H",
                "He",
                "Li",
                "Be",
                "B",
                "C",
                "N",
                "O",
                "F",
                "Ne",
                "Na",
                "Mg",
                "Al",
                "Si",
                "P",
                "S",
                "Cl",
                "Ar",
                "K",
                "Ca",
                "Sc",
                "Ti",
                "V",
                "Cr",
                "Mn",
                "Fe",
                "Co",
                "Ni",
                "Cu",
                "Zn",
                "Ga",
                "Ge",
                "As",
                "Se",
                "Br",
                "Kr",
                "Rb",
                "Sr",
                "Y",
                "Zr",
                "Nb",
                "Mo",
                "Tc",
                "Ru",
                "Rh",
                "Pd",
                "Ag",
                "Cd",
                "In",
                "Sn",
                "Sb",
                "Te",
                "I",
                "Xe",
                "Cs",
                "Ba",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Pm",
                "Sm",
                "Eu",
                "Gd",
                "Tb",
                "Dy",
                "Ho",
                "Er",
                "Tm",
                "Yb",
                "Lu",
                "Hf",
                "Ta",
                "W",
                "Re",
                "Os",
                "Ir",
                "Pt",
                "Au",
                "Hg",
                "Tl",
                "Pb",
                "Bi",
                "Po",
                "At",
                "Rn",
                "Fr",
                "Ra",
                "Ac",
                "Th",
                "Pa",
                "U",
                "Np",
                "Pu",
                "Am",
                "Cm",
                "Bk",
                "Cf",
                "Es",
                "Fm",
            ]
            return symbols[Z] if Z < len(symbols) else f"Z{Z}"

        parent_symbol = _element_symbol(self.config.Z_parent)
        daughter_symbol = _element_symbol(self.config.Z_daughter)

        header_lines = [
            f"# beta-spectrum v{__import__('beta_spectrum').__version__}",
            f"# timestamp: {datetime.now(timezone.utc).isoformat()}",
            f"# source: {source_type}",
            f"# nuclide: {parent_symbol}{self.config.A_number} -> {daughter_symbol}{self.config.A_number} (Z={self.config.Z_parent}->{self.config.Z_daughter})",
            f"# endpoint: {self.config.endpoint_MeV * 1000:.1f} keV",
            f"# transition: {self.config.transition_type}",
            f"# corrections: {', '.join(enabled)}",
            f"# e_step: {self.config.e_step_MeV:.4f} MeV",
        ]

        if self._is_multi and self.spectrum.branches:
            header_lines.append(f"# branches: {len(self.spectrum.branches)}")
            for i, branch in enumerate(self.spectrum.branches):
                header_lines.append(
                    f"#   branch_{i+1}: E₀={branch.endpoint_MeV*1000:.1f} keV, "
                    f"transition={branch.transition_type}, "
                    f"intensity={branch.intensity:.4f}"
                )
        else:
            header_lines.append("# branches: 1 (single)")

        if self.config.use_detector_response:
            header_lines.append(
                f"# detector: {self.config.detector_model} "
                f"(sigma={self.config.detector_sigma_a_keV} keV, "
                f"tail={self.config.detector_tail_fraction})"
            )
        else:
            header_lines.append("# detector: disabled")

        header_lines.append(f"# git_commit: {get_git_short_hash()}")

        # Write CSV with header
        import pandas as pd

        data: Dict[str, np.ndarray] = {
            "energy_MeV": self.energies_MeV,
            "spectrum": total,
        }

        if self._is_multi and self.spectrum.branches:
            # Multi-branch mode: total spectrum + per-branch spectra + components
            for i, (branch, branch_spec) in enumerate(
                zip(self.spectrum.branches, self.spectrum.branch_spectra)
            ):
                data[f"branch_{i+1}_spectrum"] = branch_spec(self.W)
                data[f"branch_{i+1}_intensity"] = np.full_like(
                    self.energies_MeV, branch.intensity
                )
                data[f"branch_{i+1}_endpoint"] = np.full_like(
                    self.energies_MeV, branch.endpoint_MeV
                )
                data[f"branch_{i+1}_transition"] = np.full_like(
                    self.energies_MeV, branch.transition_type, dtype=object
                )

            # Per-branch components (PhaseSpace, Radiative, etc.)
            for i, branch_spec in enumerate(self.spectrum.branch_spectra):
                branch_comps = branch_spec.calculate_components(self.W)
                for comp_name, values in branch_comps.items():
                    data[f"branch_{i+1}_{comp_name}"] = values

            # Universal components only (skip per-branch keys)
            for name, values in components.items():
                if not name.startswith("branch_"):
                    data[f"universal_{name}"] = values
        else:
            # Single-branch mode: standard export
            for name, values in components.items():
                data[name] = values

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False, float_format="%.4e")

        # Prepend header
        with open(filename, "r") as f:
            content = f.read()
        with open(filename, "w") as f:
            f.write("\n".join(header_lines) + "\n" + content)

        if self._logger:
            mode = "multi-branch" if self._is_multi else "single-branch"
            self._logger.info("CSV exported to %s (%s mode)", filename, mode)

    def get_data(self) -> Dict[str, Any]:
        """Get all numerical data for custom analysis."""
        return {
            "energies_MeV": self.energies_MeV,
            "energies_W": self.W,
            "spectrum": self.total_spectrum(normalize=True),
            "components": self.components,
            "config": self.config,
        }
