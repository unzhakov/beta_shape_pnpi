"""Tests for fitter.model."""

import numpy as np
import pytest

from beta_spectrum import BetaSpectrum, SpectrumConfig
from fitter.model import SpectrumModel


class TestSpectrumModel:
    """Test SpectrumModel."""

    def _make_model(self, with_detector: bool = False):
        """Create a test SpectrumModel."""
        config = SpectrumConfig(
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_MeV=0.294,
        )
        spectrum = BetaSpectrum.from_config(config)

        detector = None
        if with_detector:
            detector = BetaSpectrum.create_detector_from_config(
                SpectrumConfig(
                    Z_parent=43,
                    Z_daughter=44,
                    A_number=99,
                    endpoint_MeV=0.294,
                    detector_model="gaussian",
                    detector_sigma_a_keV=1.0,
                    detector_n_channels=100,
                )
            )

        model = SpectrumModel(spectrum, detector)
        return model, config

    def test_evaluate_unconvolved(self):
        model, config = self._make_model()
        energies_keV = np.linspace(0.01, 0.29, 50)
        values = model.evaluate_unconvolved(energies_keV)
        assert len(values) == 50
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0)

    def test_evaluate_no_detector(self):
        model, config = self._make_model()
        energies_keV = np.linspace(0.01, 0.29, 50)
        values = model.evaluate(energies_keV, normalization=1.0, background=0.0)
        assert len(values) == 50
        assert np.all(np.isfinite(values))

    def test_evaluate_with_normalization(self):
        model, config = self._make_model()
        energies_keV = np.linspace(0.01, 0.29, 50)
        values1 = model.evaluate(energies_keV, normalization=1.0, background=0.0)
        values2 = model.evaluate(energies_keV, normalization=2.0, background=0.0)
        np.testing.assert_allclose(values2, values1 * 2.0)

    def test_evaluate_with_background(self):
        model, config = self._make_model()
        energies_keV = np.linspace(0.01, 0.29, 50)
        values1 = model.evaluate(energies_keV, normalization=1.0, background=0.0)
        values2 = model.evaluate(energies_keV, normalization=1.0, background=5.0)
        np.testing.assert_allclose(values2, values1 + 5.0)

    def test_convolve_requires_detector(self):
        model, _ = self._make_model()
        energies_keV = np.linspace(0.01, 0.29, 50)
        with pytest.raises(RuntimeError):
            model.convolve(energies_keV)

    def test_convolve_with_detector(self):
        model, _ = self._make_model(with_detector=True)
        energies_keV = model.detector_response.channel_energies
        values = model.convolve(energies_keV)
        assert len(values) == model.detector_response.n_channels
        assert np.all(np.isfinite(values))
        # Convolution can produce small negative values due to numerical effects
        # — check that most values are non-negative
        assert np.sum(values >= 0) > len(values) * 0.9

    def test_default_params(self):
        from exp_data.spectrum import ExpSpectrum

        energies = np.linspace(0, 100, 200)
        counts = 50 * np.exp(-0.5 * ((energies - 50) / 5) ** 2) + 10
        exp_spec = ExpSpectrum(energies=energies, counts=counts)

        norm, bg = SpectrumModel.default_params(exp_spec)
        assert norm > 0
        assert bg >= 0

    def test_default_bounds(self):
        lower, upper = SpectrumModel.default_bounds(100)
        assert len(lower) == 2
        assert len(upper) == 2
        assert lower[0] < upper[0]
        assert lower[1] < upper[1]


class TestSpectrumModelKeVToW:
    """Test energy unit conversion in SpectrumModel."""

    def test_keV_to_W(self):
        model = SpectrumModel(None)  # type: ignore
        # 0 keV → W = 1.0
        W = model._keV_to_W(np.array([0.0]))
        assert np.isclose(W[0], 1.0)

        # 511 keV → W = 2.0
        W = model._keV_to_W(np.array([510.998950]))
        assert np.isclose(W[0], 2.0)

        # 294 keV → W ≈ 1.575
        W = model._keV_to_W(np.array([294.0]))
        assert np.isclose(W[0], 294.0 / 510.998950 + 1.0)
