"""Unit tests for ROOT I/O (skip if no ROOT library available)."""

from __future__ import annotations


import numpy as np
import pytest

from exp_data.spectrum import ExpSpectrum


def _has_root_io():
    """Check if any ROOT I/O library is available."""
    try:
        import ROOT  # noqa: F401

        return True
    except ImportError:
        pass
    try:
        import uproot  # noqa: F401

        return True
    except ImportError:
        return False


@pytest.fixture
def sample_spectrum():
    """Create a sample ExpSpectrum for testing."""
    return ExpSpectrum(
        energies=np.arange(100, dtype=np.float64),
        counts=np.random.poisson(10, 100).astype(np.float64),
        source="Tc99",
        run_id="test_run",
        dead_time=0.05,
        live_time=3420.0,
        date="2025-01-15",
        metadata={
            "calibration": {"type": "linear", "coefficients": [0.1, 0.4]},
            "livetime": {"fraction": 0.95},
        },
    )


@pytest.mark.skipif(
    not _has_root_io(),
    reason="Neither PyROOT nor uproot available",
)
class TestRootIOWriteRead:
    def test_round_trip(self, sample_spectrum, tmp_path):
        """Test write and read round-trip."""
        from exp_data.root_io import write_spectrum, read_spectrum

        root_path = tmp_path / "test.root"
        write_spectrum(str(root_path), sample_spectrum)
        read_back = read_spectrum(str(root_path))

        np.testing.assert_array_equal(read_back.energies, sample_spectrum.energies)
        np.testing.assert_array_equal(read_back.counts, sample_spectrum.counts)

    def test_scalar_fields(self, sample_spectrum, tmp_path):
        """Test that scalar fields survive round-trip."""
        from exp_data.root_io import write_spectrum, read_spectrum

        root_path = tmp_path / "test.root"
        write_spectrum(str(root_path), sample_spectrum)
        read_back = read_spectrum(str(root_path))

        assert read_back.dead_time == sample_spectrum.dead_time
        assert read_back.source == sample_spectrum.source
        assert read_back.run_id == sample_spectrum.run_id

    def test_metadata(self, sample_spectrum, tmp_path):
        """Test that metadata survives round-trip."""
        from exp_data.root_io import write_spectrum, read_spectrum

        root_path = tmp_path / "test.root"
        write_spectrum(str(root_path), sample_spectrum)
        read_back = read_spectrum(str(root_path))

        assert "calibration" in read_back.metadata
        assert read_back.metadata["calibration"]["type"] == "linear"
        assert read_back.metadata["livetime"]["fraction"] == 0.95
