# spectrum.py
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple
from dataclasses import dataclass

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
class SpectrumConfig:
    """
    Configuration for beta spectrum calculation.

    Supports declarative detector response specification:
    set detector_model, detector_sigma_a_keV, and related parameters
    to enable automatic detector smearing via convolve_with_detector().
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


class BetaSpectrum:
    """
    Precise beta spectrum calculator.

    Combines basic spectrum shape (fermi * phase_space) with mutplicative corrections.
    """

    def __init__(self, components: List[SpectrumComponent]):
        self.components = components
        self._component_names = [self._get_component_name(c) for c in components]

    def __call__(self, W: np.ndarray) -> np.ndarray:
        """
        Calculate total spectrum as a product of enabled components
        """
        result: np.ndarray = np.ones_like(W, dtype=float)
        for comp in self.components:
            result = result * comp(W)
        return result.astype(np.float64)

    def _get_component_name(self, comp: SpectrumComponent) -> str:
        """
        Get human-readable component names
        """
        class_name = comp.__class__.__name__
        return class_name.replace("Correction", "").replace("Function", "")

    def calculate_components(self, W: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate individual components for debugging and analysis
        """
        components_dict = {}
        for name, comp in zip(self._component_names, self.components):
            components_dict[name] = comp(W)
        return components_dict

    @classmethod
    def from_config(cls, config: SpectrumConfig) -> BetaSpectrum:
        """
        Create a BetaSpectrum from configuration
        """
        components: List[SpectrumComponent] = []
        W0 = float(T_to_W(config.endpoint_MeV))

        if config.use_phase_space:
            components.append(PhaseSpace(W0=W0, transition_type=config.transition_type))

        if config.use_fermi:
            components.append(FermiFunction(Z=config.Z_daughter, A=config.A_number))

        if config.use_finite_size:
            components.append(FiniteSizeL0(Z=config.Z_daughter, A=config.A_number))

        if config.use_charge_dist:
            components.append(
                ChargeDistributionU(Z=config.Z_daughter, A=config.A_number)
            )

        if config.use_screening:
            components.append(
                ScreeningCorrection(FermiFunction(Z=config.Z_parent, A=config.A_number))
            )

        if config.use_exchange:
            components.append(ExchangeCorrection(Z=config.Z_parent))

        if config.use_radiative:
            components.append(RadiativeCorrection(W0=W0, use_endpoint_resummation=True))

        return cls(components)

    def get_energy_grid(self, config: SpectrumConfig) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate array of energy values for spectrum calculation.
        """
        kinetic_MeV = np.arange(
            config.e_step_MeV, config.endpoint_MeV, config.e_step_MeV
        )
        W = T_to_W(kinetic_MeV)

        return W, kinetic_MeV.astype(np.float64)

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
    """

    def __init__(self, spectrum: BetaSpectrum, config: SpectrumConfig):
        """
        Initialize analyzer with a spectrum and configuration
        """
        self.spectrum = spectrum
        self.config = config
        self.W, self.energies_MeV = spectrum.get_energy_grid(config)
        self._components_cache: Optional[Dict[str, np.ndarray]] = None

    @property
    def components(self) -> Dict[str, np.ndarray]:
        if self._components_cache is None:
            self._components_cache = self.spectrum.calculate_components(self.W)
        assert self._components_cache is not None
        return self._components_cache

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

    def plot_analysis(self, save_path: Optional[str] = None) -> None:
        """
        Create visualization of the spectrum and all correction.
        """
        total = self.total_spectrum(normalize=True)
        components = self.components

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

    def export_to_csv(self, filename: str) -> None:
        total = self.total_spectrum(normalize=True)
        components = self.components

        data: Dict[str, np.ndarray] = {
            "energy_MeV": self.energies_MeV,
            "spectrum": total,
        }
        for name, values in components.items():
            data[name] = values

        import pandas as pd

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False, float_format="%.4e")
        print(f"Spectrum exported to {filename}")

    def get_data(self) -> Dict[str, Any]:
        """Get all numerical data for custom analysis."""
        return {
            "energies_MeV": self.energies_MeV,
            "energies_W": self.W,
            "spectrum": self.total_spectrum(normalize=True),
            "components": self.components,
            "config": self.config,
        }
