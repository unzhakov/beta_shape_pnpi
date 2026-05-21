"""
Physics tests for DetectorResponse — unique physical behavior.

Tests detector response shapes, resolution, and convolution physics.
"""

import numpy as np
import pytest

try:
    from beta_spectrum.components.detector_response import DetectorResponse
except ImportError:
    DetectorResponse = None  # type: ignore[misc, assignment]


@pytest.mark.skipif(DetectorResponse is None, reason="DetectorResponse not yet implemented")
class TestAnalyticalGaussianResponse:
    """Test pure Gaussian detector response."""

    def test_gaussian_peak_position(self):
        """Gaussian response should peak at the input energy."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.005)

        energy = 2.0
        response = resp.get_response(energy)
        peak_idx = np.argmax(response)
        assert np.isclose(
            channel_energies[peak_idx], energy, atol=0.01
        ), "Gaussian peak should be at input energy"

    def test_gaussian_normalization(self):
        """Gaussian response should integrate to 1 (unit area)."""
        channel_energies = np.linspace(0.0, 5.0, 8192)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.005)

        energy = 2.0
        response = resp.get_response(energy)
        integral = np.trapezoid(response, channel_energies)
        assert np.isclose(
            integral, 1.0, rtol=1e-2
        ), f"Gaussian should normalize to 1, got {integral}"

    def test_gaussian_sigma_dependence(self):
        """Larger sigma_a should produce wider Gaussian."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp_narrow = DetectorResponse(
            channel_energies, model="gaussian", sigma_a=0.002
        )
        resp_wide = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.02)

        energy = 2.0
        resp_n = resp_narrow.get_response(energy)
        resp_w = resp_wide.get_response(energy)

        half_max_n = np.max(resp_n) / 2
        half_max_w = np.max(resp_w) / 2
        fwhm_n = np.sum(resp_n > half_max_n)
        fwhm_w = np.sum(resp_w > half_max_w)
        assert fwhm_w > fwhm_n, "Larger sigma should produce wider peak"

    def test_delta_limit(self):
        """Zero sigma should produce delta function at peak channel."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(
            channel_energies, model="gaussian", sigma_a=0.0, fano_factor=0.0
        )

        energy = 2.0
        response = resp.get_response(energy)
        assert (
            np.sum(response > 0) == 1
        ), "Delta function should have single non-zero element"
        assert np.isclose(response[np.argmax(response)], 1.0)


class TestAnalyticalGaussianTailResponse:
    """Test Gaussian + exponential tail detector response."""

    def test_tail_increases_low_energy_counts(self):
        """Adding tail should produce non-zero tail contribution below peak."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp_with_tail = DetectorResponse(
            channel_energies,
            model="gaussian_tail",
            sigma_a=0.005,
            tail_fraction=0.3,
            tau=0.01,
        )

        energy = 2.0
        response = resp_with_tail.get_response(energy)

        below_peak = channel_energies < energy
        assert (
            np.sum(response[below_peak]) > 0
        ), "Response below peak should be non-zero with tail enabled"

        peak_idx = np.argmax(response)
        assert peak_idx > 0, "Peak should not be at first channel"
        assert (
            response[peak_idx] >= response[peak_idx - 1]
        ), "Response at peak should be >= response just below peak"

    def test_tail_fraction_effect(self):
        """Higher tail_fraction should produce more tail contribution."""
        channel_energies = np.linspace(0.0, 5.0, 4096)

        resp_low = DetectorResponse(
            channel_energies,
            model="gaussian_tail",
            sigma_a=0.005,
            tail_fraction=0.1,
            tau=0.01,
        )
        resp_high = DetectorResponse(
            channel_energies,
            model="gaussian_tail",
            sigma_a=0.005,
            tail_fraction=0.5,
            tau=0.01,
        )
        resp_gaussian_only = DetectorResponse(
            channel_energies, model="gaussian", sigma_a=0.005
        )

        energy = 2.0
        resp_l = resp_low.get_response(energy)
        resp_h = resp_high.get_response(energy)
        resp_g = resp_gaussian_only.get_response(energy)

        below_peak = channel_energies < energy

        tail_low = resp_l[below_peak] - (1.0 - 0.1) * resp_g[below_peak]
        tail_high = resp_h[below_peak] - (1.0 - 0.5) * resp_g[below_peak]

        assert np.sum(tail_low) > 0, "Low tail_fraction should have positive tail"
        assert np.sum(tail_high) > 0, "High tail_fraction should have positive tail"

        ratio = np.sum(tail_high) / np.sum(tail_low)
        assert ratio > 1.0, (
            f"Higher tail_fraction should have more tail contribution, "
            f"ratio={ratio:.2f}"
        )

    def test_tail_normalization(self):
        """Gaussian+tail response should still integrate to 1."""
        channel_energies = np.linspace(0.0, 5.0, 8192)
        resp = DetectorResponse(
            channel_energies,
            model="gaussian_tail",
            sigma_a=0.005,
            tail_fraction=0.2,
            tau=0.01,
        )

        energy = 2.0
        response = resp.get_response(energy)
        integral = np.trapezoid(response, channel_energies)
        assert np.isclose(
            integral, 1.0, rtol=1e-2
        ), f"Gaussian+tail should normalize to 1, got {integral}"


class TestResolutionSigma:
    """Test energy-dependent resolution sigma(E)."""

    def test_sigma_increases_with_energy(self):
        """σ(E) should increase with energy for sigma_b > 0."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(
            channel_energies, model="gaussian", sigma_a=0.0, sigma_b=0.01
        )

        sigma_low = resp._resolution_sigma(0.5)
        sigma_high = resp._resolution_sigma(4.0)
        assert sigma_high > sigma_low, "σ should increase with √E"

    def test_sigma_a_only(self):
        """With sigma_b=0, σ should be constant."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(
            channel_energies,
            model="gaussian",
            sigma_a=0.01,
            sigma_b=0.0,
            fano_factor=0.0,
        )

        sigma_low = resp._resolution_sigma(0.5)
        sigma_high = resp._resolution_sigma(4.0)
        assert np.isclose(sigma_low, sigma_high), "σ should be constant when sigma_b=0"

    def test_fano_contribution(self):
        """Fano factor should add energy-dependent contribution to σ."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp_no_fano = DetectorResponse(
            channel_energies,
            model="gaussian",
            sigma_a=0.0,
            sigma_b=0.0,
            fano_factor=0.0,
        )
        resp_fano = DetectorResponse(
            channel_energies,
            model="gaussian",
            sigma_a=0.0,
            sigma_b=0.0,
            fano_factor=0.12,
        )

        sigma_low = resp_no_fano._resolution_sigma(0.5)
        sigma_fano = resp_fano._resolution_sigma(0.5)
        assert sigma_fano > sigma_low, "Fano factor should increase σ"


class TestTabulatedResponse:
    """Test tabulated (MC-simulated) detector response."""

    def test_tabulated_interpolation(self):
        """Tabulated response should interpolate between calibration energies."""
        channel_energies = np.linspace(0.0, 5.0, 256)
        calib_energies = np.array([1.0, 2.0, 3.0, 4.0])

        n_channels = len(channel_energies)
        response_matrix = np.zeros((len(calib_energies), n_channels))
        for i, e in enumerate(calib_energies):
            idx = np.argmin(np.abs(channel_energies - e))
            response_matrix[i, idx] = 1.0

        resp = DetectorResponse(
            channel_energies=channel_energies,
            response_matrix=response_matrix,
            calibration_energies=calib_energies,
        )

        response = resp.get_response(2.5)
        assert np.any(response > 0), "Interpolated response should be non-zero"

    def test_tabulated_outside_range(self):
        """Response outside calibration range should be zero."""
        channel_energies = np.linspace(0.0, 5.0, 256)
        calib_energies = np.array([1.0, 2.0, 3.0])

        response_matrix = np.zeros((3, len(channel_energies)))
        for i, e in enumerate(calib_energies):
            idx = np.argmin(np.abs(channel_energies - e))
            response_matrix[i, idx] = 1.0

        resp = DetectorResponse(
            channel_energies=channel_energies,
            response_matrix=response_matrix,
            calibration_energies=calib_energies,
        )

        response = resp.get_response(0.5)
        assert np.all(response == 0), "Response below calibration range should be zero"

        response = resp.get_response(4.0)
        assert np.all(response == 0), "Response above calibration range should be zero"


class TestConvolution:
    """Test spectrum convolution with detector response."""

    def test_convolution_preserves_total_counts(self):
        """Total counts should be approximately preserved after convolution."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.005)

        W = np.linspace(0.0, 5.0, 1000)
        spectrum = np.zeros_like(W)
        idx = np.argmin(np.abs(W - 2.5))
        spectrum[idx] = 1.0 / (W[1] - W[0])

        convolved = resp.convolve(W, spectrum, normalize=True)
        total_theoretical = np.trapezoid(spectrum, W)
        total_convolved = np.trapezoid(convolved, channel_energies)
        assert np.isclose(
            total_convolved, total_theoretical, rtol=0.05
        ), f"Total counts: theory={total_theoretical:.4f}, convolved={total_convolved:.4f}"

    def test_convolution_widens_peak(self):
        """Convolution should widen a narrow peak."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.02)

        W = np.linspace(0.0, 5.0, 2000)
        spectrum = np.zeros_like(W)
        idx = np.argmin(np.abs(W - 2.5))
        spectrum[idx] = 1.0 / (W[1] - W[0])

        convolved = resp.convolve(W, spectrum, normalize=True)
        peak_width = np.sum(convolved > np.max(convolved) / 2)
        assert peak_width > 1, "Convolved peak should have finite width"

    def test_convolution_beta_spectrum(self):
        """Convolution of a beta spectrum should produce a smooth curve."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.01)

        W0 = 5.0
        W = np.linspace(1.0, W0 - 0.01, 500)
        p = np.sqrt(W**2 - 1)
        spectrum = p * W * (W0 - W) ** 2

        convolved = resp.convolve(W, spectrum, normalize=True)
        assert len(convolved) == len(channel_energies)
        assert np.all(convolved >= 0), "Convolved spectrum should be non-negative"
        assert np.max(convolved) > 0, "Convolved spectrum should have non-zero peak"

    def test_convolve_batch_same_as_convolve(self):
        """Batch convolution should match standard convolution for analytical models."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(
            channel_energies,
            model="gaussian_tail",
            sigma_a=0.01,
            tail_fraction=0.1,
            tau=0.01,
        )

        W = np.linspace(1.0, 4.9, 300)
        p = np.sqrt(W**2 - 1)
        spectrum = p * W * (4.9 - W) ** 2

        result_standard = resp.convolve(W, spectrum, normalize=True)
        result_batch = resp.convolve_batch(W, spectrum, normalize=True)

        assert np.allclose(
            result_standard, result_batch, rtol=1e-5
        ), "Batch and standard convolution should match for analytical models"
