# tests/test_detector_response.py — Detector response unit tests

"""
Why test this?
--------------
The detector response function is the bridge between theoretical beta spectra
and measured experimental data. It models how a mono-energetic electron beam
appears in the detector's energy spectrum — typically a Gaussian peak with
possible low-energy tailing due to charge collection inefficiencies.

We test:

 1. **Analytical response shapes**: Gaussian core, Gaussian+tail, normalization
 2. **Tabulated response**: interpolation from MC-calculated response matrix
 3. **Convolution**: theoretical spectrum → predicted measured spectrum
 4. **Resolution model**: energy-dependent σ(E) behavior
 5. **Edge cases**: zero resolution, zero tail fraction, boundary energies
"""

import numpy as np

from beta_spectrum.components.detector_response import DetectorResponse


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

        # FWHM: count points where response > half max
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

        # Below peak: the tail contribution should be non-zero
        below_peak = channel_energies < energy
        # The response below peak comes from both the Gaussian tail (left side)
        # and the exponential tail. Both should be non-zero.
        assert (
            np.sum(response[below_peak]) > 0
        ), "Response below peak should be non-zero with tail enabled"

        # The response at the peak should be at or above nearby energies
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

        # The tail contribution is the difference between gaussian_tail
        # and pure gaussian (scaled by tail_fraction).
        # Higher tail_fraction should have larger tail contribution.
        # Tail contribution = f_tail * T(E), where T is the tail function.
        # Since T is only non-zero for E < E0, the tail contribution
        # below peak is proportional to f_tail.
        tail_low = resp_l[below_peak] - (1.0 - 0.1) * resp_g[below_peak]
        tail_high = resp_h[below_peak] - (1.0 - 0.5) * resp_g[below_peak]

        # Both should have positive tail contribution
        assert np.sum(tail_low) > 0, "Low tail_fraction should have positive tail"
        assert np.sum(tail_high) > 0, "High tail_fraction should have positive tail"

        # Higher tail_fraction should have proportionally more tail contribution
        # (tail_high / tail_low should be roughly 0.5/0.1 = 5)
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

        # Simple diagonal response matrix: each calib energy peaks at its channel
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

        # Response at 2.5 should interpolate between 2.0 and 3.0 rows
        response = resp.get_response(2.5)
        assert np.any(response > 0), "Interpolated response should be non-zero"

    def test_tabulated_validation_wrong_columns(self):
        """Should reject response matrix with wrong column count."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        calib_energies = np.array([1.0, 2.0])
        response_matrix = np.zeros((2, 200))  # 200 columns, but 100 channels

        try:
            DetectorResponse(
                channel_energies=channel_energies,
                response_matrix=response_matrix,
                calibration_energies=calib_energies,
            )
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

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

        # Energy below calibration range
        response = resp.get_response(0.5)
        assert np.all(response == 0), "Response below calibration range should be zero"

        # Energy above calibration range
        response = resp.get_response(4.0)
        assert np.all(response == 0), "Response above calibration range should be zero"


class TestConvolution:
    """Test spectrum convolution with detector response."""

    def test_convolution_preserves_total_counts(self):
        """Total counts should be approximately preserved after convolution."""
        channel_energies = np.linspace(0.0, 5.0, 4096)
        resp = DetectorResponse(channel_energies, model="gaussian", sigma_a=0.005)

        # Mono-energetic line at W=2.5
        W = np.linspace(0.0, 5.0, 1000)
        spectrum = np.zeros_like(W)
        idx = np.argmin(np.abs(W - 2.5))
        spectrum[idx] = 1.0 / (W[1] - W[0])  # Delta function approximation

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

        # Very narrow peak
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

        # Simple phase space spectrum ∝ p·W·(W0-W)²
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
        spectrum = np.zeros(50)  # Wrong length
        try:
            resp.convolve(W, spectrum)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass

    def test_non_1d_input(self):
        """Should reject non-1D inputs."""
        channel_energies = np.linspace(0.0, 5.0, 100)
        resp = DetectorResponse(channel_energies, model="gaussian")
        W = np.array([[0.0, 5.0]])  # 2D
        spectrum = np.array([[1.0, 1.0]])
        try:
            resp.convolve(W, spectrum)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass


class TestBetaSpectrumIntegration:
    """Integration tests between BetaSpectrum and DetectorResponse."""

    def test_beta_spectrum_convolve_with_detector(self):
        """BetaSpectrum should convolve with detector response."""
        from beta_spectrum import BetaSpectrum, SpectrumConfig, DetectorResponse
        from beta_spectrum.utils import T_to_W

        config = SpectrumConfig(
            Z_parent=43, Z_daughter=44, A_number=99, endpoint_MeV=0.294
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        # Detector channel energies must be in m_e units (total energy)
        W0 = T_to_W(0.294)
        detector = DetectorResponse.from_gaussian_params(
            channel_energy_range=(1.0, W0 + 0.05),
            n_channels=512,
            sigma_a=0.003,
            sigma_b=0.0,
            tail_fraction=0.0,
            model="gaussian",
            fano_factor=0.0,
        )

        convolved = spectrum.convolve_with_detector(detector, W=W, config=config)

        assert convolved.shape == (512,)
        assert np.all(convolved >= 0), "Convolved spectrum should be non-negative"
        assert np.sum(convolved) > 0, "Convolved spectrum should have non-zero area"

    def test_convolved_spectrum_has_wider_peak(self):
        """Convolved spectrum should be wider than theoretical spectrum."""
        from beta_spectrum import BetaSpectrum, SpectrumConfig, DetectorResponse
        from beta_spectrum.utils import T_to_W

        config = SpectrumConfig(
            Z_parent=43, Z_daughter=44, A_number=99, endpoint_MeV=0.294
        )
        spectrum = BetaSpectrum.from_config(config)
        W, _ = spectrum.get_energy_grid(config)

        # Theoretical spectrum
        theoretical = spectrum(W)

        # Detector response with finite resolution
        W0 = T_to_W(0.294)
        detector = DetectorResponse.from_gaussian_params(
            channel_energy_range=(1.0, W0 + 0.05),
            n_channels=512,
            sigma_a=0.005,
            sigma_b=0.0,
            tail_fraction=0.0,
            model="gaussian",
            fano_factor=0.0,
        )

        convolved = spectrum.convolve_with_detector(detector, W=W, config=config)

        # The convolved spectrum should be smoother/wider
        theoretical_width = np.sum(theoretical > np.max(theoretical) * 0.5)
        convolved_width = np.sum(convolved > np.max(convolved) * 0.5)
        assert (
            convolved_width >= theoretical_width
        ), "Convolved spectrum should be at least as wide as theoretical"

    def test_analyzer_convolved_spectrum(self):
        """BetaSpectrumAnalyzer should compute convolved spectrum."""
        from beta_spectrum import (
            BetaSpectrum,
            SpectrumConfig,
            DetectorResponse,
            BetaSpectrumAnalyzer,
        )
        from beta_spectrum.utils import T_to_W

        config = SpectrumConfig(
            Z_parent=43, Z_daughter=44, A_number=99, endpoint_MeV=0.294
        )
        spectrum = BetaSpectrum.from_config(config)
        analyzer = BetaSpectrumAnalyzer(spectrum, config)

        W0 = T_to_W(0.294)
        detector = DetectorResponse.from_gaussian_params(
            channel_energy_range=(1.0, W0 + 0.05),
            n_channels=256,
            sigma_a=0.003,
            sigma_b=0.0,
            tail_fraction=0.0,
            model="gaussian",
            fano_factor=0.0,
        )

        convolved = analyzer.convolved_spectrum(detector, normalize=True)

        assert convolved.shape == (256,)
        assert np.isclose(
            np.trapezoid(convolved, detector.channel_energies), 1.0, rtol=0.01
        ), "Normalized convolved spectrum should integrate to 1"
