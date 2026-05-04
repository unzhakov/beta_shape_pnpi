"""Tests for logging infrastructure."""

import logging
import os

import pytest

from beta_spectrum.logging_utils import (
    LoggingConfig,
    _find_git_root,
    get_git_short_hash,
    get_logger,
    setup_logging,
)


class TestLoggingConfig:
    """Test LoggingConfig dataclass defaults and behavior."""

    def test_default_level(self):
        config = LoggingConfig()
        assert config.level == "WARNING"

    def test_default_log_file(self):
        config = LoggingConfig()
        assert config.log_file is None

    def test_custom_level(self):
        config = LoggingConfig(level="DEBUG")
        assert config.level == "DEBUG"

    def test_custom_log_file(self):
        config = LoggingConfig(log_file="test.log")
        assert config.log_file == "test.log"


class TestSetupLogging:
    """Test logging setup and configuration."""

    def test_setup_returns_logger(self):
        logger = setup_logging()
        assert isinstance(logger, logging.Logger)
        assert logger.name == "beta_spectrum"

    def test_setup_sets_level(self):
        config = LoggingConfig(level="DEBUG")
        logger = setup_logging(config)
        assert logger.level == logging.DEBUG

    def test_setup_info_level(self):
        config = LoggingConfig(level="INFO")
        logger = setup_logging(config)
        assert logger.level == logging.INFO

    def test_setup_creates_console_handler(self):
        config = LoggingConfig()
        logger = setup_logging(config)
        handlers = [h for h in logger.handlers if isinstance(h, logging.StreamHandler)]
        assert len(handlers) == 1

    def test_setup_creates_file_handler(self, tmp_path):
        log_file = str(tmp_path / "test.log")
        config = LoggingConfig(log_file=log_file)
        logger = setup_logging(config)
        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert os.path.exists(log_file)

    def test_setup_prevents_duplicate_handlers(self):
        config = LoggingConfig()
        logger1 = setup_logging(config)
        logger2 = setup_logging(config)
        assert logger1 is logger2
        assert len(logger1.handlers) == 1  # Only console handler

    def test_setup_with_file_prevents_duplicates(self, tmp_path):
        log_file = str(tmp_path / "test2.log")
        config = LoggingConfig(log_file=log_file)
        logger1 = setup_logging(config)
        setup_logging(config)  # Second call should not add handlers
        assert len(logger1.handlers) == 2  # Console + File


class TestGetLogger:
    """Test getting named loggers."""

    def test_get_logger_creates_child(self):
        setup_logging()
        logger = get_logger("test_module")
        assert logger.name == "beta_spectrum.test_module"
        assert logger.parent is not None

    def test_get_logger_nested(self):
        setup_logging()
        logger = get_logger("spectrum.components")
        assert logger.name == "beta_spectrum.spectrum.components"


class TestGitShortHash:
    """Test git hash retrieval."""

    def test_returns_short_hash(self):
        result = get_git_short_hash()
        assert len(result) == 7
        assert all(c in "0123456789abcdef" for c in result)

    def test_returns_custom_length(self):
        result = get_git_short_hash(length=5)
        assert len(result) == 5

    def test_returns_unknown_outside_repo(self):
        # Patch _find_git_root to return a non-existent path
        import unittest.mock

        with unittest.mock.patch(
            "beta_spectrum.logging_utils._find_git_root",
            return_value="/nonexistent/path",
        ):
            result = get_git_short_hash()
            assert result == "unknown"

    def test_find_git_root(self):
        root = _find_git_root()
        assert os.path.exists(os.path.join(root, ".git"))


class TestCSVMetadataHeader:
    """Test CSV export metadata header."""

    def test_header_contains_version(self):
        from beta_spectrum.spectrum import (
            BetaSpectrum,
            BetaSpectrumAnalyzer,
            SpectrumConfig,
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F2",
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            filename = f.name

        analyzer.export_to_csv(filename, source_type="paceENSDF")

        with open(filename, "r") as f:
            lines = f.readlines()

        assert any("# beta-spectrum v" in line for line in lines)
        assert any("0.3.0" in line for line in lines)
        os.unlink(filename)

    def test_header_contains_timestamp(self):
        from beta_spectrum.spectrum import (
            BetaSpectrum,
            BetaSpectrumAnalyzer,
            SpectrumConfig,
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F2",
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            filename = f.name

        analyzer.export_to_csv(filename, source_type="json")

        with open(filename, "r") as f:
            content = f.read()

        assert "# timestamp:" in content
        os.unlink(filename)

    def test_header_contains_nuclide_info(self):
        from beta_spectrum.spectrum import (
            BetaSpectrum,
            BetaSpectrumAnalyzer,
            SpectrumConfig,
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F2",
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            filename = f.name

        analyzer.export_to_csv(filename, source_type="paceENSDF")

        with open(filename, "r") as f:
            content = f.read()

        assert "# nuclide: 43->44, A=99" in content
        os.unlink(filename)

    def test_header_contains_git_commit(self):
        from beta_spectrum.spectrum import (
            BetaSpectrum,
            BetaSpectrumAnalyzer,
            SpectrumConfig,
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F2",
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            filename = f.name

        analyzer.export_to_csv(filename, source_type="cli")

        with open(filename, "r") as f:
            content = f.read()

        assert "# git_commit:" in content
        os.unlink(filename)

    def test_header_contains_corrections(self):
        from beta_spectrum.spectrum import (
            BetaSpectrum,
            BetaSpectrumAnalyzer,
            SpectrumConfig,
        )

        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
            transition_type="F2",
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            filename = f.name

        analyzer.export_to_csv(filename, source_type="paceENSDF")

        with open(filename, "r") as f:
            content = f.read()

        assert "# corrections:" in content
        assert "phase_space" in content
        assert "fermi" in content
        os.unlink(filename)


class TestCLILogging:
    """Test CLI logging integration."""

    def test_verbose_flag_sets_info(self):
        """Test that -v sets INFO level."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--nuclide", "Tc99", "--dry-run"])
        assert args.verbose == 0

    def test_double_verbose_flag(self):
        """Test that -vv sets DEBUG level."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--nuclide", "Tc99", "--dry-run", "-vv"])
        assert args.verbose == 2

    def test_quiet_flag(self, capsys):
        """Test that -q suppresses output."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--nuclide", "Tc99", "--dry-run", "-q"])
        assert args.quiet is True

    def test_log_file_flag(self, capsys):
        """Test that --log-file sets log file path."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        args = parser.parse_args(
            ["--nuclide", "Tc99", "--dry-run", "--log-file", "test.log"]
        )
        assert args.log_file == "test.log"

    def test_dry_run_flag(self, capsys):
        """Test that --dry-run is recognized."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        args = parser.parse_args(["--nuclide", "Tc99", "--dry-run"])
        assert args.dry_run is True

    def test_version_flag(self, capsys):
        """Test that --version prints version."""
        from beta_spectrum.cli import _build_parser

        parser = _build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--version"])
        assert exc_info.value.code == 0
