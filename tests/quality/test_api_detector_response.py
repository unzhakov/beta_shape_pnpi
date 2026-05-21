"""
API / quality tests for DetectorResponse.

Tests factory methods, repr, error handling.
"""

import numpy as np
import pytest

try:
    from beta_spectrum.components.detector_response import DetectorResponse
except ImportError:
    DetectorResponse = None  # type: ignore[misc, assignment]


@pytest.mark.skipif(DetectorResponse is None, reason="DetectorResponse not yet implemented")
class TestFactoryMethods:
    """Test convenience factory methods."""

    def test_from_gaussian_params(self):
        """from_gaussian_params should create correct DetectorResponse."""
        resp = DetectorResponse.from_gaussian_params(
            channel_energy_range=(0.0, 5.0),
            n_channels=1024,
            sigma_a=0.005,
            sigma_b=0.01,
            tail_fraction=0.2,
            tau=0.01,
            model="gaussian_tail",
            fano_factor=0.12,
        )

        assert resp.n_channels == 1024
        assert np.isclose(resp.channel_energies[0], 0.0)
        assert np.isclose(resp.channel_energies[-1], 5.0)
        assert resp._mode == "analytical"
        assert resp.tail_fraction == 0.2

    def test_from_mc_simulation(self):
        """from_mc_simulation should create tabulated DetectorResponse."""
        channel_energies = np.linspace(0.0, 5.0, 512)
        calib_energies = np.array([1.0, 2.0, 3.0])
        response_matrix = np.zeros((3, 512))
        for i, e in enumerate(calib_energies):
            idx = np.argmin(np.abs(channel_energies - e))
            response_matrix[i, idx] = 1.0

        resp = DetectorResponse.from_mc_simulation(
            channel_energies=channel_energies,
            response_matrix=response_matrix,
            calibration_energies=calib_energies,
        )

        assert resp._mode == "tabulated"
        assert resp.n_channels == 512
        assert len(resp.calibration_energies) == 3


@pytest.mark.skipif(DetectorResponse is None, reason="DetectorResponse not yet implemented")
class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_repr_analytical(self):
        """__repr__ should show analytical mode info."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.01)
        repr_str = repr(resp)
        assert "analytical" in repr_str
        assert "gaussian" in repr_str

    def test_repr_tabulated(self):
        """__repr__ should show tabulated mode info."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        calib_energies = np.array([1.0, 2.0])
        response_matrix = np.zeros((2, 100))
        resp = DetectorResponse(
            channel_energies=channel_energies,
            response_matrix=response_matrix,
            calibration_energies=calib_energies,
        )
        repr_str = repr(resp)
        assert "tabulated" in repr_str

    def test_invalid_model(self):
        """Should reject unknown model names."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(channel_energies, model="unknown_model")
        try:
            resp.get_response(2.0)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Unknown model" in str(e)

    def test_convolve_mismatched_lengths(self):
        """Should reject mismatched W and spectrum_values lengths."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(channel_energies, model="gaussian")
        W = np.linspace(0.0, 5.0, 100)
        spectrum = np.zeros(50)
        try:
            resp.convolve(W, spectrum)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

    def test_non_1d_input(self):
        """Should reject non-1D inputs."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(channel_energies, model="gaussian")
        W = np.array([[0.0, 5.0]])
        spectrum = np.array([[1.0, 1.0]])
        try:
            resp.convolve(W, spectrum)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

    def test_tabulated_validation_wrong_columns(self):
        """Should reject response matrix with wrong column count."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        calib_energies = np.array([1.0, 2.0])
        response_matrix = np.zeros((2, 200))

        try:
            DetectorResponse(
                channel_energies=channel_energies,
                response_matrix=response_matrix,
                calibration_energies=calib_energies,
            )
            assert False, "Should have raised ValueError"
        except ValueError:
            pass
