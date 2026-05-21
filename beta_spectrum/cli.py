"""
Command-line interface for beta-spectrum calculations.

Usage examples:
    # From paceENSDF nuclear data
    bs_pnpi --nuclide Tc99 --output spectrum.csv

    # From custom JSON input
    bs_pnpi --input custom_input.json --output spectrum.csv

    # Full analysis with plotting
    bs_pnpi --nuclide Co60 --plot analysis.png --output spectrum.csv

    # Verbose with log file
    bs_pnpi --nuclide Tc99 --output spectrum.csv -vv --log-file /tmp/calc.log

    # Dry run to inspect resolved configuration
    bs_pnpi --nuclide Tc99 --dry-run
"""

from __future__ import annotations

import argparse
import sys

from beta_spectrum.logging_utils import LoggingConfig, setup_logging
from beta_spectrum.spectrum import SpectrumConfig


def _build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog="bs_pnpi",
        description="Beta spectrum calculation toolkit — "
        "calculate high-precision beta decay energy spectra.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  bs_pnpi --nuclide Tc99 --output spectrum.csv\n"
            "  bs_pnpi --input custom_input.json --output spectrum.csv\n"
            "  bs_pnpi --nuclide Tc99 --plot spectrum.png --output spectrum.csv\n"
            "  bs_pnpi --nuclide Tc99 --detector gaussian --output convolved.csv\n"
            "  bs_pnpi --nuclide Tc99 --output spectrum.csv -vv --log-file /tmp/calc.log\n"
            "  bs_pnpi --nuclide Tc99 --dry-run\n"
        ),
    )

    # Version
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__import__('beta_spectrum').__version__}",
    )

    # Source selection (mutually exclusive)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--nuclide",
        type=str,
        help="Nuclide symbol for paceENSDF lookup (e.g. 'Tc99', 'Co60'). "
        "Retrieves decay parameters from the ENSDF database via paceENSDF.",
    )
    source_group.add_argument(
        "--input",
        type=str,
        metavar="FILE",
        help="Path to a JSON input file with custom parameters.",
    )

    # Common parameters
    parser.add_argument(
        "--e-step",
        type=float,
        default=0.001,
        help="Energy step size in MeV (default: 0.001). All energies are in MeV throughout the toolkit.",
    )

    # Output
    parser.add_argument(
        "--output",
        type=str,
        metavar="FILE",
        help="Output CSV file path.",
    )

    # Plotting
    parser.add_argument(
        "--plot",
        type=str,
        metavar="FILE",
        help="Save analysis plot to the specified path.",
    )

    # Logging and output control
    log_group = parser.add_argument_group("logging and output")
    log_group.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity: -v=INFO (shows source, components, and workflow), "
        "-vv=DEBUG (detailed internals of all spectrum components). Repeatable.",
    )
    log_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all terminal output (only errors shown).",
    )
    log_group.add_argument(
        "--log-file",
        type=str,
        metavar="PATH",
        help="Write log messages to the specified file path (absolute or relative). "
        "Log file is always opened in append mode at DEBUG level.",
    )
    log_group.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate input and display resolved configuration without computing spectrum.",
    )
    log_group.add_argument(
        "--intensity-cutoff",
        type=float,
        default=0.0,
        metavar="FRAC",
        help="Minimum branch intensity fraction to include in the spectrum "
        "(default: 0.0 = include all branches). "
        "Branches below this threshold are excluded from the calculation and plots.",
    )

    return parser


def _run(args: argparse.Namespace) -> None:
    """Execute the CLI command."""
    import time

    from beta_spectrum.logging_utils import get_logger
    from beta_spectrum.nuclear_data import (
        create_config_from_source,
        get_decay_info_from_paceENSDF,
        json_to_config,
        load_json_input,
    )
    from beta_spectrum.spectrum import (
        BetaSpectrum,
        BranchConfig,
        SpectrumConfig,
    )
    from beta_spectrum.visualize import BetaSpectrumAnalyzer

    # Determine logging level
    if args.quiet:
        log_level = "WARNING"
    elif args.verbose >= 2:
        log_level = "DEBUG"
    elif args.verbose >= 1:
        log_level = "INFO"
    else:
        log_level = "WARNING"

    # Setup logging
    log_config = LoggingConfig(
        level=log_level,
        log_file=args.log_file,
    )
    logger = setup_logging(log_config)
    logger = get_logger("cli")

    logger.info("bs_pnpi starting...")
    logger.debug("Verbosity: %s | Log file: %s", log_level, args.log_file or "none")

    # Build config
    source_type: str
    config: SpectrumConfig
    if args.nuclide:
        logger.info("Source: paceENSDF nuclide=%s", args.nuclide)
        decay_info = get_decay_info_from_paceENSDF(args.nuclide, "beta_minus")

        # Normalize branch intensities to sum to 1.0
        if decay_info.branches:
            total_intensity = sum(b.intensity for b in decay_info.branches)
            if total_intensity > 0:
                branches = [
                    BranchConfig(
                        endpoint_MeV=(
                            (decay_info.endpoint_MeV - b.level_energy_keV / 1000.0)
                            if b.level_energy_keV > 0
                            else decay_info.endpoint_MeV
                        ),
                        transition_type=b.transition_type,
                        intensity=b.intensity / total_intensity,
                    )
                    for b in decay_info.branches
                    if b.intensity / total_intensity >= args.intensity_cutoff
                ]
                config = SpectrumConfig(
                    Z_parent=decay_info.Z_parent,
                    Z_daughter=decay_info.Z_daughter,
                    A_number=decay_info.A_number,
                    endpoint_MeV=decay_info.endpoint_MeV,
                    transition_type=decay_info.transition_type,
                    e_step_MeV=args.e_step,
                    branches=branches,
                    intensity_cutoff=args.intensity_cutoff,
                )
            else:
                config = create_config_from_source(
                    "paceENSDF",
                    nuclide=args.nuclide,
                    e_step_MeV=args.e_step,
                )
        else:
            config = create_config_from_source(
                "paceENSDF",
                nuclide=args.nuclide,
                e_step_MeV=args.e_step,
            )
        source_type = "paceENSDF"
    elif args.input:
        logger.info("Source: JSON file=%s", args.input)
        data = load_json_input(args.input)
        config = json_to_config(data)
        source_type = "json"
    else:
        logger.error("No source specified (use --nuclide or --input)")
        sys.exit(1)

    # Dry run: display resolved config and exit
    if args.dry_run:
        _print_dry_run_output(config, source_type)
        return

    # Create and calculate spectrum
    start_time = time.time()
    mode_str = ""
    if config.branches and len(config.branches) > 1:
        mode_str = f" (multi-branch: {len(config.branches)} branches)"
    logger.info(
        "Creating BetaSpectrum for %s, endpoint=%.3f MeV%s",
        _decay_notation(config.Z_parent, config.Z_daughter, config.A_number),
        config.endpoint_MeV,
        mode_str,
    )

    spectrum = BetaSpectrum.from_config(config, logger=logger)
    _W, _kinetic_MeV = spectrum.get_energy_grid(config)
    _ = spectrum(_W)

    elapsed = time.time() - start_time
    logger.info(
        "Spectrum calculation done in %.3fs (%d points)", elapsed, len(_kinetic_MeV)
    )

    # Export to CSV
    if args.output:
        analyzer = BetaSpectrumAnalyzer(spectrum, config, logger=logger)
        analyzer.export_to_csv(args.output, source_type=source_type)
        logger.info("Spectrum exported to %s", args.output)

    # Generate plot
    if args.plot:
        analyzer = BetaSpectrumAnalyzer(spectrum, config, logger=logger)
        analyzer.plot_analysis(
            save_path=args.plot,
            show_components=(args.verbose >= 2),
        )
        logger.info("Analysis plot saved to %s", args.plot)

    logger.info("Done.")


def _print_dry_run_output(config, source_type: str) -> None:
    """Print resolved configuration for dry-run mode."""
    from beta_spectrum.logging_utils import get_git_short_hash

    enabled = []
    if config.use_phase_space:
        enabled.append("phase_space")
    if config.use_fermi:
        enabled.append("fermi")
    if config.use_screening:
        enabled.append("screening")
    if config.use_finite_size:
        enabled.append("finite_size")
    if config.use_charge_dist:
        enabled.append("charge_dist")
    if config.use_radiative:
        enabled.append("radiative")
    if config.use_exchange:
        enabled.append("exchange")

    print("=== Dry Run: Resolved Configuration ===")
    print(f"Package:    beta-spectrum v{__import__('beta_spectrum').__version__}")
    print(f"Git commit: {get_git_short_hash()}")
    print(f"Source:     {source_type}")
    print(
        f"Nuclide:    {_decay_notation(config.Z_parent, config.Z_daughter, config.A_number)}"
    )
    print(f"Endpoint:   {config.endpoint_MeV * 1000:.1f} keV")
    print(f"Transition: {config.transition_type}")
    print(f"e_step:     {config.e_step_MeV:.4f} MeV")
    print(f"Corrections: {', '.join(enabled) if enabled else 'none'}")

    if config.branches and len(config.branches) > 1:
        print(f"\nMulti-branch mode: {len(config.branches)} branches")
        if config.intensity_cutoff > 0:
            print(
                f"Intensity cutoff: {config.intensity_cutoff:.4f} ({config.intensity_cutoff*100:.1f}%)"
            )
        for i, branch in enumerate(config.branches):
            print(
                f"  Branch {i+1}: E₀={branch.endpoint_MeV*1000:.1f} keV, "
                f"transition={branch.transition_type}, "
                f"intensity={branch.intensity:.4f}"
            )
    else:
        print("Mode: single-branch")

    if config.use_detector_response:
        print(
            f"Detector:   {config.detector_model} (sigma={config.detector_sigma_a_keV} keV, "
            f"tail={config.detector_tail_fraction})"
        )
    else:
        print("Detector:   disabled")
    print("========================================")


def _element_symbol(Z: int) -> str:
    """Get element symbol from atomic number."""
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


def _nuclide_symbol(Z: int, A: int) -> str:
    """Get nuclide notation like 'Tc99' from atomic and mass numbers."""
    return f"{_element_symbol(Z)}{A}"


def _decay_notation(Z_parent: int, Z_daughter: int, A: int) -> str:
    """Get decay notation like 'Tc99 -> Ru99' from atomic and mass numbers."""
    return f"{_nuclide_symbol(Z_parent, A)} -> {_nuclide_symbol(Z_daughter, A)}"


def main() -> None:
    """Main entry point for the CLI."""
    parser = _build_parser()
    args = parser.parse_args()

    try:
        _run(args)
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
