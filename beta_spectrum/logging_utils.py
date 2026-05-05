"""
Logging infrastructure for beta-spectrum package.

Provides a centralized logging setup with:
- Configurable verbosity levels (WARNING, INFO, DEBUG)
- Optional file output
- Consistent formatting across all modules
- Opt-in design: logging is disabled by default for library use,
  enabled via setup_logging() for CLI/workflow use

Usage (library):
    # No logging by default — silent, zero overhead

Usage (CLI):
    from beta_spectrum.logging_utils import setup_logging, get_logger
    setup_logging(level="INFO", log_file="calc.log")
    logger = get_logger("spectrum")
    logger.info("Starting calculation...")
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class LoggingConfig:
    """Configuration for logging setup."""

    level: str = "WARNING"
    """Logging level: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'."""

    log_file: Optional[str] = None
    """Path to log file. If None, only stdout is used."""

    log_format: str = "[%(asctime)s] [%(levelname)-8s] [%(name)s] %(message)s"
    """Log message format string."""

    date_format: str = "%H:%M:%S"
    """Date/time format for timestamps."""

    _handlers_configured: bool = field(default=False, repr=False)
    """Internal flag to prevent double configuration."""


def setup_logging(config: Optional[LoggingConfig] = None) -> logging.Logger:
    """
    Configure the root logger for the beta_spectrum package.

    This is a one-shot setup — calling it multiple times with different
    configs will not reconfigure existing handlers. Call setup_logging()
    once at the start of your application.

    Parameters
    ----------
    config : LoggingConfig, optional
        Logging configuration. Uses defaults if None.

    Returns
    -------
    logging.Logger
        The root logger for the beta_spectrum package.
    """
    if config is None:
        config = LoggingConfig()

    logger = logging.getLogger("beta_spectrum")

    # Prevent double configuration
    if config._handlers_configured:
        return logger

    logger.setLevel(getattr(logging, config.level.upper(), logging.WARNING))

    # Remove any existing handlers to avoid duplicates on re-init
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, config.level.upper(), logging.WARNING))
    console_handler.setFormatter(
        logging.Formatter(config.log_format, config.date_format)
    )
    logger.addHandler(console_handler)

    # File handler (optional)
    if config.log_file:
        file_handler = logging.FileHandler(config.log_file, mode="a")
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_handler.setFormatter(
            logging.Formatter(config.log_format, config.date_format)
        )
        logger.addHandler(file_handler)

    config._handlers_configured = True
    return logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a named logger for the beta_spectrum package.

    Parameters
    ----------
    name : str
        Logger name (e.g., "spectrum", "fitter", "nuclear_data").

    Returns
    -------
    logging.Logger
        A logger named "beta_spectrum.{name}".

    Examples
    --------
    >>> logger = get_logger("spectrum")
    >>> logger.info("Creating spectrum components")
    """
    return logging.getLogger(f"beta_spectrum.{name}")


def get_git_short_hash(length: int = 7) -> str:
    """
    Get the short git commit hash.

    Parameters
    ----------
    length : int
        Number of characters to return (default: 7).

    Returns
    -------
    str
        Short git hash, or "unknown" if not in a git repo.
    """
    import subprocess

    try:
        result = subprocess.run(
            ["git", "rev-parse", f"--short={length}", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=_find_git_root(),
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _find_git_root() -> str:
    """Find the root directory of the git repository."""
    from pathlib import Path

    current = Path(__file__).resolve()
    while current != current.parent:
        if (current / ".git").exists():
            return str(current)
        current = current.parent
    return str(current)
