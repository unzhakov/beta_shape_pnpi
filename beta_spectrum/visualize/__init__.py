"""
Visualize module — BetaSpectrumAnalyzer class.

Visualization and plotting tools for BetaSpectrum.
All nuclides (single- and multi-branch) are treated uniformly as one or more
branches — single-branch is just a special case with one branch.
"""

from __future__ import annotations

import logging
import pandas as pd
from datetime import datetime, timezone
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt

from beta_spectrum.logging_utils import get_git_short_hash
from beta_spectrum.utils import T_to_W

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from beta_spectrum.spectrum import BetaSpectrum, SpectrumConfig
    from beta_spectrum.components.detector_response import DetectorResponse


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

        # Unified branch view: single-branch treated as one-branch multi-branch
        if not self._is_multi and self.spectrum.branch_spectra == []:
            self._unified_branches: List[_UnifiedBranch] = [
                _UnifiedBranch(
                    endpoint_MeV=config.endpoint_MeV,
                    transition_type=config.transition_type,
                    intensity=1.0,
                )
            ]
        else:
            self._unified_branches = [
                _UnifiedBranch(
                    endpoint_MeV=b.endpoint_MeV,
                    transition_type=b.transition_type,
                    intensity=b.intensity,
                )
                for b in self.spectrum.branches
            ]

    @property
    def components(self) -> Dict[str, np.ndarray]:
        if self._components_cache is None:
            raw = self.spectrum.calculate_components(self.W)
            # Normalize keys: single-branch has bare names,
            # multi-branch has branch_0., branch_1., etc.
            # We keep the raw keys for export_to_csv compatibility.
            self._components_cache = raw
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
        """Calculate the total spectrum."""
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

    # ---------------------------------------------------------------------------
    # Unified branch view
    # ---------------------------------------------------------------------------

    @property
    def _n_branches(self) -> int:
        """Number of branches (unified view)."""
        return len(self._unified_branches)

    def _get_branch_spectrum(self, i: int) -> np.ndarray:
        """Get the spectrum array for branch i."""
        if i < len(self.spectrum.branch_spectra):
            return self.spectrum.branch_spectra[i]
        # Single-branch: the spectrum itself is the only branch
        return self.spectrum(self.W)

    def _get_component_key(self, i: int, name: str) -> str:
        """Get the canonical key for a per-branch component."""
        if self._is_multi:
            return f"branch_{i}.{name}"
        return name

    def _has_component(self, i: int, name: str) -> bool:
        """Check if branch i has a given component."""
        key = self._get_component_key(i, name)
        return key in self.components

    def _get_component(self, i: int, name: str) -> np.ndarray:
        """Get a per-branch component array."""
        key = self._get_component_key(i, name)
        return self.components[key]

    # ---------------------------------------------------------------------------
    # Plotting helpers
    # ---------------------------------------------------------------------------

    def _build_figure_header(self, fig: "Figure") -> None:
        """Add a compact header with nuclide info, e_step, commit and timestamp outside axes."""
        parent = self._element_symbol(self.config.Z_parent)
        daughter = self._element_symbol(self.config.Z_daughter)
        e_step_keV = self.config.e_step_MeV * 1000
        header = (
            f"{parent}{self.config.A_number} -> {daughter}{self.config.A_number}  "
            f"|  E₀={self.config.endpoint_MeV * 1000:.1f} keV  "
            f"|  ΔE={e_step_keV:.1f} keV"
        )
        fig.text(
            0.01, 0.98, header,
            transform=fig.transFigure, fontsize=9,
            verticalalignment="top", horizontalalignment="left",
            fontfamily="monospace",
        )
        id_text = f"{get_git_short_hash(6)}  |  {self._timestamp}"
        fig.text(
            0.01, 0.01, id_text,
            transform=fig.transFigure, fontsize=8,
            verticalalignment="bottom", horizontalalignment="left",
            fontfamily="monospace", color="gray",
        )

    def _build_figure_title(self, fig: "Figure") -> str:
        """Build a unified figure title."""
        parent = self._element_symbol(self.config.Z_parent)
        daughter = self._element_symbol(self.config.Z_daughter)
        n = self.config.A_number
        title = f"Beta-decay: {n}{parent} -> {n}{daughter}"
        if self._n_branches > 1:
            title += f" ({self._n_branches} branches)"
        return title

    def _format_branch_label(self, i: int) -> str:
        """Abbreviated branch label for legend."""
        branch = self._unified_branches[i]
        e0 = branch.endpoint_MeV * 1000
        intensity_pct = branch.intensity * 100
        return f"Br. {i+1}: {e0:.1f} keV ({intensity_pct:.1f}%)"

    # ---------------------------------------------------------------------------
    # Plotting — unified multi-branch approach
    # ---------------------------------------------------------------------------

    def plot_analysis(
        self,
        save_path: Optional[str] = None,
        show_components: bool = True,
    ) -> None:
        """
        Create visualization of the spectrum and all correction factors.

        Parameters
        ----------
        save_path : str, optional
            Path to save the figure.
        show_components : bool
            If True, show full debug view.
            If False, show only the spectrum plot.
        """
        total = self.total_spectrum(normalize=True)
        components = self.components
        self._timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        if show_components:
            self._plot_debug_view(total, components, save_path)
        else:
            self._plot_spectra_view(total, components, save_path)

    def _plot_spectra_view(
        self,
        total: np.ndarray,
        components: Dict[str, np.ndarray],
        save_path: Optional[str],
    ) -> None:
        """Plot total spectrum with individual branch decomposition."""
        fig, ax = plt.subplots(figsize=(12, 7))

        branch_colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf",
        ]

        e_step_keV = self.config.e_step_MeV * 1000

        # Compute branch contributions (universal * intensity * branch_spectrum)
        branch_contribs = []
        for i in range(self._n_branches):
            universal_product = np.ones_like(self.W)
            for name in [
                "Fermi", "Screening", "FiniteSizeL0",
                "ChargeDistributionU", "Exchange",
            ]:
                key = self._get_component_key(i, name)
                if key in components:
                    universal_product *= components[key]
            bs = self._get_branch_spectrum(i)
            contrib = universal_product * self._unified_branches[i].intensity * bs
            branch_contribs.append(contrib)

        total_raw = sum(branch_contribs)

        ax.plot(
            self.energies_MeV, total_raw, "k-", lw=2.5,
            label="Total spectrum", zorder=10,
        )

        for i, contrib in enumerate(branch_contribs):
            color = branch_colors[i % len(branch_colors)]
            ax.plot(
                self.energies_MeV, contrib, color=color, lw=1.5,
                alpha=0.7, label=self._format_branch_label(i),
            )

        ax.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=12)
        ax.set_ylabel(
            f"Normalized Counts / {e_step_keV:.1f} keV", fontsize=12
        )
        fig.suptitle(
            self._build_figure_title(fig), fontsize=14, fontweight="bold"
        )
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend(
            loc="best", fontsize=9,
            ncol=min(2, self._n_branches + 1),
        )

        self._build_figure_header(fig)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

    def _plot_debug_view(
        self,
        total: np.ndarray,
        components: Dict[str, np.ndarray],
        save_path: Optional[str],
    ) -> None:
        """Debug view: vertical layout with 4 panels."""
        n_branches = self._n_branches
        e_step_keV = self.config.e_step_MeV * 1000

        fig, axes = plt.subplots(4, 1, figsize=(12, 16), sharex=True)
        fig.suptitle(
            self._build_figure_title(fig), fontsize=14, fontweight="bold"
        )

        branch_colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf",
        ]

        # Panel 1: Total spectrum + per-branch spectra
        branch_contribs = []
        for i in range(n_branches):
            universal_product = np.ones_like(self.W)
            for name in [
                "Fermi", "Screening", "FiniteSizeL0",
                "ChargeDistributionU", "Exchange",
            ]:
                key = self._get_component_key(i, name)
                if key in components:
                    universal_product *= components[key]
            bs = self._get_branch_spectrum(i)
            contrib = universal_product * self._unified_branches[i].intensity * bs
            branch_contribs.append(contrib)

        total_raw = sum(branch_contribs)

        ax = axes[0]
        ax.plot(
            self.energies_MeV, total_raw, "k-", lw=2,
            label="Total spectrum", zorder=10,
        )
        for i, contrib in enumerate(branch_contribs):
            color = branch_colors[i % len(branch_colors)]
            ax.plot(
                self.energies_MeV, contrib, color=color, lw=1.5,
                alpha=0.7, label=self._format_branch_label(i),
            )
        ax.set_yscale("log")
        ax.set_ylabel(
            f"Normalized Counts / {e_step_keV:.1f} keV", fontsize=10
        )
        ax.set_title("Spectrum", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best", ncol=min(2, n_branches + 1))

        # Panel 2: Fermi function + all phase spaces
        ax = axes[1]
        if "Fermi" in components:
            fermi_vals = components["Fermi"]
            fermi_scaled = fermi_vals / np.max(fermi_vals)
            ax.plot(
                self.energies_MeV, fermi_scaled, "r-", lw=1.5,
                label="Fermi", alpha=0.8,
            )
        for i in range(n_branches):
            color = branch_colors[i % len(branch_colors)]
            branch = self._unified_branches[i]
            key = self._get_component_key(i, "PhaseSpace")
            if key in components:
                ps_vals = components[key]
                W0_branch = T_to_W(branch.endpoint_MeV)
                mask = self.W <= W0_branch
                ps_masked = np.where(mask, ps_vals, 0.0)
                ps_max = np.max(ps_masked)
                if ps_max > 0:
                    ps_scaled = ps_masked / ps_max
                else:
                    ps_scaled = ps_masked
                ax.plot(
                    self.energies_MeV, ps_scaled, color=color, lw=1.5,
                    alpha=0.7, label=f"PS Br. {i+1}",
                )
        ax.set_ylabel("Normalized Factor", fontsize=10)
        ax.set_title("Fermi & Phase Space", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")

        # Panel 3: Universal W0-independent corrections
        ax = axes[2]
        universal_components = [
            "Screening", "FiniteSizeL0",
            "ChargeDistributionU", "Exchange",
        ]
        colors = ["green", "orange", "purple", "brown"]
        for name, color in zip(universal_components, colors):
            if name in components:
                values = components[name]
                values_scaled = values / np.max(values)
                ax.plot(
                    self.energies_MeV, values_scaled, color=color, lw=1.5,
                    label=name, alpha=0.8,
                )
        ax.set_ylabel("Normalized Factor", fontsize=10)
        ax.set_title("Universal Corrections", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")

        # Panel 4: W0-dependent corrections (Radiative for each branch)
        ax = axes[3]
        for i in range(n_branches):
            color = branch_colors[i % len(branch_colors)]
            key = self._get_component_key(i, "Radiative")
            if key in components:
                rad_vals = components[key]
                ax.plot(
                    self.energies_MeV, rad_vals, color=color, lw=1.5,
                    alpha=0.7, label=f"Br. {i+1}",
                )
        ax.set_ylabel("Radiative Correction", fontsize=10)
        ax.set_xlabel(r"Electron kinetic energy $E$ [MeV]", fontsize=10)
        ax.set_title("Radiative Corrections", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best", ncol=min(2, n_branches + 1))

        self._build_figure_header(fig)

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Figure saved to {save_path}")

    def _element_symbol(self, Z: int) -> str:
        """Return chemical symbol for atomic number Z."""
        symbols = [
            "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
            "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
            "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
            "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
            "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
            "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
            "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
            "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
            "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
            "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
        ]
        return symbols[Z] if Z < len(symbols) else f"Z{Z}"

    # ---------------------------------------------------------------------------
    # Export
    # ---------------------------------------------------------------------------

    def export_to_csv(self, filename: str, source_type: str = "unknown") -> None:
        """Export spectrum data to CSV with metadata header."""
        total = self.total_spectrum(normalize=True)
        components = self.components

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

        parent_symbol = self._element_symbol(self.config.Z_parent)
        daughter_symbol = self._element_symbol(self.config.Z_daughter)

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

        # Branch info — unified view
        for i, branch in enumerate(self._unified_branches):
            header_lines.append(
                f"#   branch_{i+1}: E₀={branch.endpoint_MeV*1000:.1f} keV, "
                f"transition={branch.transition_type}, "
                f"intensity={branch.intensity:.4f}"
            )

        if self.config.use_detector_response:
            header_lines.append(
                f"# detector: {self.config.detector_model} "
                f"(sigma={self.config.detector_sigma_a_keV} keV, "
                f"tail={self.config.detector_tail_fraction})"
            )
        else:
            header_lines.append("# detector: disabled")

        header_lines.append(f"# git_commit: {get_git_short_hash()}")

        data: Dict[str, np.ndarray] = {
            "energy_MeV": self.energies_MeV,
            "spectrum": total,
        }

        # Per-branch spectra and components
        for i in range(self._n_branches):
            data[f"branch_{i+1}_spectrum"] = self._get_branch_spectrum(i)
            for comp_name in ["PhaseSpace", "Radiative"]:
                key = self._get_component_key(i, comp_name)
                if key in components:
                    data[f"branch_{i+1}_{comp_name}"] = components[key]

        # Universal components (no prefix)
        for name, values in components.items():
            if not name.startswith("branch_"):
                data[name] = values

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False, float_format="%.4e")

        with open(filename, "r") as f:
            content = f.read()
        with open(filename, "w") as f:
            f.write("\n".join(header_lines) + "\n" + content)

        if self._logger:
            self._logger.info("CSV exported to %s", filename)

    def get_data(self) -> Dict[str, Any]:
        """Get all numerical data for custom analysis."""
        return {
            "energies_MeV": self.energies_MeV,
            "energies_W": self.W,
            "spectrum": self.total_spectrum(normalize=True),
            "components": self.components,
            "config": self.config,
        }


@dataclass
class _UnifiedBranch:
    """Unified branch view for both single- and multi-branch modes."""
    endpoint_MeV: float
    transition_type: str
    intensity: float
