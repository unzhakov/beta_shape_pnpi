"""Tests for exp_data.calibration."""

import numpy as np
import pytest

from exp_data.calibration import CalibrationResult, EnergyCalibrator


class TestCalibrationResult:
    """Test CalibrationResult."""

    def test_linear_call(self):
        # polyfit returns [slope, intercept] = [0.5, 5.0]
        # polyval([0.5, 5.0], ch) = 0.5*ch + 5.0
        result = CalibrationResult(
            coefficients=np.array([0.5, 5.0]),  # E = 0.5*ch + 5.0
            order=1,
            chi2_per_dof=0.01,
            n_points=2,
            covariance=np.eye(2),
        )
        np.testing.assert_allclose(result(0), 5.0)
        np.testing.assert_allclose(result(10), 10.0)

    def test_quadratic_call(self):
        # polyfit returns [a2, a1, a0] = [0.01, 0.1, 1.0]
        # polyval([0.01, 0.1, 1.0], ch) = 0.01*ch² + 0.1*ch + 1.0
        result = CalibrationResult(
            coefficients=np.array([0.01, 0.1, 1.0]),
            order=2,
            chi2_per_dof=0.01,
            n_points=3,
            covariance=np.eye(3),
        )
        np.testing.assert_allclose(result(0), 1.0, rtol=1e-10)
        # ch=10: 0.01*100 + 0.1*10 + 1.0 = 1 + 1 + 1 = 3
        np.testing.assert_allclose(result(10), 3.0)

    def test_invert(self):
        result = CalibrationResult(
            coefficients=np.array([0.5, 5.0]),  # E = 0.5*ch + 5.0
            order=1,
            chi2_per_dof=0.01,
            n_points=2,
            covariance=np.eye(2),
        )
        inv = result.invert()
        # E=10 => ch = (10-5)/0.5 = 10
        np.testing.assert_allclose(inv(10.0), 10.0, rtol=0.01)

    def test_resolution_at(self):
        result = CalibrationResult(
            coefficients=np.array([1.0, 0.0]),  # E = 1.0*ch + 0
            order=1,
            chi2_per_dof=0.01,
            n_points=2,
            covariance=np.eye(2),
        )
        # dE/dch = 1.0, resolution = 1.0 * 1.0 = 1.0 keV
        res = result.resolution_at(50.0, sigma_channel=1.0)
        np.testing.assert_allclose(res, 1.0)

    def test_resolution_at_quadratic(self):
        result = CalibrationResult(
            coefficients=np.array([0.02, 1.0, 0.0]),  # E = 0.02*ch² + 1.0*ch
            order=2,
            chi2_per_dof=0.01,
            n_points=3,
            covariance=np.eye(3),
        )
        # dE/dch = 2*0.02*ch + 1.0 = 0.04*ch + 1.0
        # at ch=50: dE/dch = 0.04*50 + 1.0 = 3.0
        expected = 0.04 * 50.0 + 1.0
        res = result.resolution_at(50.0, sigma_channel=1.0)
        np.testing.assert_allclose(res, expected)


class TestEnergyCalibrator:
    """Test EnergyCalibrator."""

    def test_not_calibrated_initially(self):
        cal = EnergyCalibrator()
        assert not cal.is_calibrated

    def test_get_linear_calibration(self):
        cal = EnergyCalibrator()
        result = cal.get_linear_calibration(
            low_energy_keV=0.0,
            high_energy_keV=100.0,
            low_channel=0,
            high_channel=100,
        )
        assert cal.is_calibrated
        assert result.order == 1
        np.testing.assert_allclose(result(50), 50.0)

    def test_apply_calibration_with_points(self):
        cal = EnergyCalibrator()
        cal.add_calibration_point(0, 0.0)
        cal.add_calibration_point(50, 50.0)
        cal.add_calibration_point(100, 100.0)

        result = cal.fit_calibration(order=1)
        assert result.order == 1
        # Perfect fit on collinear points
        np.testing.assert_allclose(result(50), 50.0, rtol=1e-6)

    def test_channel_to_energy_raises_without_calibration(self):
        cal = EnergyCalibrator()
        with pytest.raises(RuntimeError):
            cal.channel_to_energy(50)

    def test_energy_to_channel_raises_without_calibration(self):
        cal = EnergyCalibrator()
        with pytest.raises(RuntimeError):
            cal.energy_to_channel(50.0)

    def test_full_channel_mapping(self):
        """Calibrate using full channel-to-energy mapping."""
        # Simulate: E = 0.5 + 0.1*ch + 0.001*ch²
        n_channels = 200
        channels = np.arange(n_channels, dtype=np.float64)
        energies = 0.5 + 0.1 * channels + 0.001 * channels**2

        cal = EnergyCalibrator()
        result = cal.fit_calibration(order=2, channel_energies=energies)

        assert result.order == 2
        np.testing.assert_allclose(result(channels), energies, rtol=1e-10)

    def test_quadratic_calibration(self):
        cal = EnergyCalibrator()
        cal.add_calibration_point(0, 0.0)
        cal.add_calibration_point(50, 52.0)
        cal.add_calibration_point(100, 108.0)  # slight quadratic deviation

        result = cal.fit_calibration(order=2)
        assert result.order == 2
