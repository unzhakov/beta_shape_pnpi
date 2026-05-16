# spectrum.py
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple
from dataclasses import dataclass

import numpy as np

from beta_spectrum.base import SpectrumComponent
from beta_spectrum.components.phase_space import PhaseSpace
from beta_spectrum.components.fermi import FermiFunction
from beta_spectrum.components.finite_size import FiniteSizeL0, ChargeDistributionU
from beta_spectrum.components.screening import ScreeningCorrection
from beta_spectrum.components.exchange import ExchangeCorrection
from beta_spectrum.components.radiative import RadiativeCorrection
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


