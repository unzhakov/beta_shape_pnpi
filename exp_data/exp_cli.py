"""
Command-line interface for experimental data processing.

Usage examples:
    # Process a single spectrum
    bs_exp process --raw A1.DAT --calibration CAM1.DAT --output output.root

    # Process a run sequence
    bs_exp process --run-dir ./data/raw/ --prefix A --calibration CAM1.DAT

    # Determine live-time from a pulser calibration
    bs_exp livetime --file CAM1.DAT --frequency 100.0 --real-time 3600.0

    # Export calibrated spectrum to ROOT
    bs_exp export --input spectrum.root --output processed.root

    # Dry run to inspect pipeline configuration
    bs_exp process --raw A1.DAT --calibration CAM1.DAT --dry-run
"""

from __future__ import annotations

import argparse
import sys


def _build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        prog="bs_exp",
        description="Experimental data processing toolkit - "
        "calibrate, process, and export beta decay spectra.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  bs_exp process --raw A1.DAT --calibration CAM1.DAT --output output.root\n"
            "  bs_exp process --run-dir ./data/ --prefix A --calibration CAM1.DAT\n"
            "  bs_exp livetime --file CAM1.DAT --frequency 100.0 --real-time 3600.0\n"
            "  bs_exp export --input spectrum.root --output processed.root\n"
            "  bs_exp process --raw A1.DAT --calibration CAM1.DAT --dry-run\n"
        ),
    )

    # Version
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.4.0",
    )

    # Subcommands
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # --- process subcommand ---
    process_parser = subparsers.add_parser(
        "process",
        help="Process raw binary spectra through the full pipeline",
        description=(
            "Read raw binary spectra, subtract cumulative runs, "
            "apply energy calibration, determine live-time, and export."
        ),
    )
    process_parser.add_argument(
        "--raw",
        type=str,
        metavar="FILE",
        help="Path to the raw binary spectrum file (A*.DAT).",
    )
    process_parser.add_argument(
        "--cumulative",
        type=str,
        metavar="FILE",
        help="Path to the previous cumulative spectrum for subtraction.",
    )
    process_parser.add_argument(
        "--run-dir",
        type=str,
        metavar="DIR",
        help="Directory containing cumulative run files (for sequence processing).",
    )
    process_parser.add_argument(
        "--prefix",
        type=str,
        default="A",
        help="Filename prefix for run files (default: 'A').",
    )
    process_parser.add_argument(
        "--calibration",
        type=str,
        required=True,
        metavar="FILE",
        help="Path to the Am-241 calibration binary file (CAM1.DAT).",
    )
    process_parser.add_argument(
        "--pulser-frequency",
        type=float,
        default=None,
        metavar="HZ",
        help="Pulser injection frequency in Hz (for live-time determination).",
    )
    process_parser.add_argument(
        "--real-time",
        type=float,
        default=None,
        metavar="SEC",
        help="Real acquisition time in seconds (for live-time determination).",
    )
    process_parser.add_argument(
        "--output",
        type=str,
        metavar="FILE",
        help="Output ROOT file path.",
    )
    process_parser.add_argument(
        "--run-id",
        type=str,
        default=None,
        help="Run identifier prefix.",
    )
    process_parser.add_argument(
        "--source",
        type=str,
        default=None,
        help="Source identifier (e.g., 'Tc99').",
    )

    # --- livetime subcommand ---
    livetime_parser = subparsers.add_parser(
        "livetime",
        help="Determine live-time fraction from a pulser calibration spectrum",
        description=(
            "Read a pulser calibration spectrum, identify the pulser "
            "peak, and compute the live-time fraction."
        ),
    )
    livetime_parser.add_argument(
        "--file",
        type=str,
        required=True,
        metavar="FILE",
        help="Path to the pulser calibration binary file.",
    )
    livetime_parser.add_argument(
        "--frequency",
        type=float,
        required=True,
        metavar="HZ",
        help="Pulser injection frequency in Hz.",
    )
    livetime_parser.add_argument(
        "--real-time",
        type=float,
        required=True,
        metavar="SEC",
        help="Real acquisition time in seconds.",
    )
    livetime_parser.add_argument(
        "--output",
        type=str,
        metavar="FILE",
        help="Output ROOT file with livetime metadata.",
    )

    # --- export subcommand ---
    export_parser = subparsers.add_parser(
        "export",
        help="Export calibrated ExpSpectrum to ROOT format",
        description="Export an ExpSpectrum (from ROOT or dict) to ROOT file format.",
    )
    export_parser.add_argument(
        "--input",
        type=str,
        required=True,
        metavar="FILE",
        help="Input file (ROOT or JSON dict).",
    )
    export_parser.add_argument(
        "--output",
        type=str,
        required=True,
        metavar="FILE",
        help="Output ROOT file path.",
    )

    # --- dry-run (global, added to all subparsers) ---
    for sub in [process_parser, livetime_parser, export_parser]:
        sub.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate input and display configuration without processing.",
        )

    # Logging and output (global)
    log_group = parser.add_argument_group("logging and output")
    log_group.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity: -v=INFO, -vv=DEBUG. Repeatable.",
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
        help="Write log messages to the specified file.",
    )

    return parser


def _run(args: argparse.Namespace) -> None:
    """Execute the CLI command."""
    # Determine logging level
    if args.quiet:
        log_level = "WARNING"
    elif args.verbose >= 2:
        log_level = "DEBUG"
    elif args.verbose >= 1:
        log_level = "INFO"
    else:
        log_level = "WARNING"

    # Setup logging (minimal - reuse beta_spectrum.logging_utils if available)
    try:
        from beta_spectrum.logging_utils import setup_logging, get_logger

        log_config = type(
            "LoggingConfig", (), {"level": log_level, "log_file": args.log_file}
        )()
        logger = setup_logging(log_config)
        logger = get_logger("exp")
    except ImportError:
        import logging

        logging.basicConfig(level=getattr(logging, log_level))
        logger = logging.getLogger("exp")

    logger.info("bs_exp starting...")

    if args.command == "process":
        _cmd_process(args, logger)
    elif args.command == "livetime":
        _cmd_livetime(args, logger)
    elif args.command == "export":
        _cmd_export(args, logger)
    else:
        parser = _build_parser()
        parser.print_help()
        sys.exit(1)


def _cmd_process(args, logger) -> None:
    """Process subcommand handler."""
    from exp_data.pipeline import process_spectrum, process_run_sequence
    from exp_data.root_io import write_spectrum

    # Check PyROOT/uproot availability
    from exp_data.root_io import _PYROOT_AVAILABLE, _UPTOOL_AVAILABLE

    if not _PYROOT_AVAILABLE and not _UPTOOL_AVAILABLE:
        logger.error("No ROOT I/O library available. Install: pip install uproot")
        sys.exit(1)

    if args.dry_run:
        logger.info("Dry run: pipeline configuration validated")
        logger.info("  Calibration: %s", args.calibration)
        logger.info("  Pulser frequency: %s Hz", args.pulser_frequency_hz or "N/A")
        logger.info("  Real time: %s s", args.real_time or "N/A")
        return

    if args.raw:
        # Single spectrum
        result = process_spectrum(
            raw_path=args.raw,
            calibration_path=args.calibration,
            cumulative_path=args.cumulative,
            run_id=args.run_id,
            source=args.source,
            pulser_frequency_hz=args.pulser_frequency_hz,
            real_time_sec=args.real_time,
        )

        logger.info("Pipeline complete:")
        logger.info(
            "  Validation: %s", "PASSED" if result.validation_passed else "FAILED"
        )
        for msg in result.validation_messages:
            logger.info("  %s", msg)
        logger.info("  Calibration chi2/dof: %.4f", result.calibration.chi2_per_dof)

        if result.livetime:
            logger.info(
                "  Live-time fraction: %.4f", result.livetime.live_time_fraction
            )

        if args.output:
            write_spectrum(args.output, result.spectrum)
            logger.info("  Exported to %s", args.output)

    elif args.run_dir:
        # Run sequence
        results = process_run_sequence(
            run_dir=args.run_dir,
            calibration_path=args.calibration,
            prefix=args.prefix,
            run_id=args.run_id,
            source=args.source,
            pulser_frequency_hz=args.pulser_frequency_hz,
            real_time_sec=args.real_time,
        )

        logger.info("Processed %d runs", len(results))
        for i, result in enumerate(results):
            logger.info(
                "  Run %d: %s (chi2/dof=%.4f)",
                i + 1,
                "PASSED" if result.validation_passed else "FAILED",
                result.calibration.chi2_per_dof,
            )

        if args.output:
            # Export last result
            write_spectrum(args.output, results[-1].spectrum)
            logger.info("  Exported last run to %s", args.output)
    else:
        logger.error("Specify --raw FILE or --run-dir DIR")
        sys.exit(1)


def _cmd_livetime(args, logger) -> None:
    """Livetime subcommand handler."""
    from exp_data.livetime import LivetimeDetermination
    from exp_data.root_io import write_spectrum
    from exp_data.spectrum import ExpSpectrum

    result = LivetimeDetermination.from_file(
        filepath=args.file,
        pulser_frequency_hz=args.frequency,
        real_time_sec=args.real_time,
        source="pulser",
    )

    logger.info("Live-time determination:")
    logger.info("  Live-time fraction: %.4f", result.live_time_fraction)
    logger.info("  Observed pulser rate: %.1f Hz", result.observed_pulser_rate)
    logger.info(
        "  Pulser peak: %.2f keV (sigma=%.2f keV, chi2/dof=%.2f)",
        result.pulser_peak_center_keV,
        result.pulser_peak_sigma_keV,
        result.pulser_peak_chi2,
    )

    if args.output:
        # Create a dummy spectrum with livetime metadata
        import numpy as np

        dummy = ExpSpectrum(
            energies=np.array([0.0]),
            counts=np.array([0.0]),
            errors=np.array([0.0]),
            metadata={
                "livetime": {
                    "fraction": result.live_time_fraction,
                    "pulser_frequency_hz": result.pulser_frequency_hz,
                    "observed_pulser_rate": result.observed_pulser_rate,
                }
            },
            live_time=result.real_time_sec * result.live_time_fraction,
            source="pulser",
        )
        write_spectrum(args.output, dummy)
        logger.info("  Livetime metadata exported to %s", args.output)


def _cmd_export(args, logger) -> None:
    """Export subcommand handler."""
    from exp_data.root_io import write_spectrum, read_spectrum

    # Read input spectrum
    if str(args.input).endswith(".root"):
        spectrum = read_spectrum(args.input)
    elif str(args.input).endswith(".json"):
        import json
        from exp_data.spectrum import ExpSpectrum

        with open(args.input) as f:
            data = json.load(f)
        spectrum = ExpSpectrum.from_dict(data)
    else:
        logger.error("Unsupported input format. Use .root or .json")
        sys.exit(1)

    # Write output
    write_spectrum(args.output, spectrum)
    logger.info("Exported to %s", args.output)


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
