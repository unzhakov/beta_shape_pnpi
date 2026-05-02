"""
Command-line interface for beta-spectrum calculations.

Usage examples:
    # From paceENSDF nuclear data
    bs_pnpi --nuclide Tc99 --output spectrum.csv

    # From custom JSON input
    bs_pnpi --input custom_input.json --output spectrum.csv

    # Full analysis with plotting
    bs_pnpi --nuclide Co60 --mode BM --plot analysis.png --output spectrum.csv

    # With detector response
    bs_pnpi --nuclide Tc99 --detector sigma=1.0 --output convolved.csv
"""

from __future__ import annotations

import argparse
import sys
from typing import Any, Dict


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
            "  bs_pnpi --nuclide Tc99 --detector sigma=1.0 --output convolved.csv\n"
        ),
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
        "--mode",
        type=str,
        default="beta_minus",
        choices=["beta_minus", "beta_plus", "ec"],
        help="Decay mode for paceENSDF source (default: beta_minus).",
    )
    parser.add_argument(
        "--decay-index",
        type=int,
        default=None,
        help="Decay index for paceENSDF (0 = ground state, default: auto-select).",
    )
    parser.add_argument(
        "--transition-type",
        type=str,
        default=None,
        choices=["A", "F1", "F1U", "F2", "F2U", "F3", "F3U", "F4"],
        help="Override transition type (forbiddenness classification).",
    )
    parser.add_argument(
        "--e-step",
        type=float,
        default=0.001,
        help="Energy step size in MeV (default: 0.001).",
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

    # Detector response
    detector_group = parser.add_argument_group("detector response")
    detector_group.add_argument(
        "--detector",
        type=str,
        default=None,
        metavar="MODEL",
        help="Enable detector response convolution. "
        "Specify model: 'gaussian', 'gaussian_tail', 'tikhonov', 'tabulated'.",
    )
    detector_group.add_argument(
        "--sigma",
        type=float,
        default=1.0,
        metavar="KEV",
        help="Resolution parameter sigma in keV (default: 1.0).",
    )
    detector_group.add_argument(
        "--tau",
        type=float,
        default=5.0,
        metavar="KEV",
        help="Tail decay constant tau in keV (default: 5.0).",
    )
    detector_group.add_argument(
        "--tail-fraction",
        type=float,
        default=0.0,
        metavar="FRAC",
        help="Low-energy tail fraction (0.0–1.0, default: 0.0).",
    )

    return parser


def _parse_detector_args(args: argparse.Namespace) -> Dict[str, Any]:
    """Parse detector-related arguments into a config dict."""
    if not args.detector:
        return {}

    return {
        "use_detector_response": True,
        "detector_model": args.detector,
        "detector_sigma_a_keV": args.sigma,
        "detector_tau_keV": args.tau,
        "detector_tail_fraction": args.tail_fraction,
    }


def _run(args: argparse.Namespace) -> None:
    """Execute the CLI command."""
    from beta_spectrum.nuclear_data import (
        create_config_from_source,
        json_to_config,
        load_json_input,
    )
    from beta_spectrum.spectrum import BetaSpectrum, BetaSpectrumAnalyzer

    # Build config
    detector_kwargs = _parse_detector_args(args)

    if args.nuclide:
        config = create_config_from_source(
            "paceENSDF",
            nuclide=args.nuclide,
            decay_type=args.mode,
            decay_index=args.decay_index,
            e_step_MeV=args.e_step,
            use_detector_response=detector_kwargs.get("use_detector_response", False),
            detector_model=detector_kwargs.get("detector_model", "gaussian"),
            detector_sigma_a_keV=detector_kwargs.get("detector_sigma_a_keV", 1.0),
            detector_tau_keV=detector_kwargs.get("detector_tau_keV", 5.0),
            detector_tail_fraction=detector_kwargs.get("detector_tail_fraction", 0.0),
        )
    elif args.input:
        data = load_json_input(args.input)
        if args.transition_type:
            data["transition_type"] = args.transition_type
        config = json_to_config(data)

    # Override transition type if explicitly provided
    if args.transition_type:
        config.transition_type = args.transition_type

    # Create and calculate spectrum
    spectrum = BetaSpectrum.from_config(config)
    _W, _kinetic_MeV = spectrum.get_energy_grid(config)
    _ = spectrum(_W)

    # Export to CSV
    if args.output:
        analyzer = BetaSpectrumAnalyzer(spectrum, config)
        analyzer.export_to_csv(args.output)
        print(f"Spectrum exported to {args.output}")

    # Generate plot
    if args.plot:
        analyzer = BetaSpectrumAnalyzer(spectrum, config)
        analyzer.plot_analysis(save_path=args.plot)
        print(f"Analysis plot saved to {args.plot}")


def main() -> None:
    """Main entry point for the CLI."""
    parser = _build_parser()
    args = parser.parse_args()

    try:
        _run(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
