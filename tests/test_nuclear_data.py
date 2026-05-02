"""Tests for nuclear_data module (paceENSDF integration and JSON input)."""

import json
import tempfile
from pathlib import Path

import pytest

# ruff: noqa: E402, ANN201, ANN001
from beta_spectrum.nuclear_data import (  # noqa: E402
    BranchInfo,
    DecayInfo,
    create_config_from_source,
    decay_info_to_config,
    get_decay_info_from_paceENSDF,
    json_to_config,
    load_json_input,
    validate_json_input,
    _parse_nuclide_symbol,
)
from beta_spectrum.spectrum import SpectrumConfig


def _has_paceENSDF() -> bool:  # noqa: ANN201
    """Check if paceENSDF is available."""
    try:
        import paceENSDF  # noqa: F401

        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# _parse_nuclide_symbol
# ---------------------------------------------------------------------------


class TestParseNuclideSymbol:
    """Test nuclide symbol parsing."""

    def test_tc99(self):
        element, Z, A = _parse_nuclide_symbol("Tc99")
        assert element == "Tc"
        assert Z == 43
        assert A == 99

    def test_co60(self):
        element, Z, A = _parse_nuclide_symbol("Co60")
        assert element == "Co"
        assert Z == 27
        assert A == 60

    def test_he3(self):
        element, Z, A = _parse_nuclide_symbol("He3")
        assert element == "He"
        assert Z == 2
        assert A == 3

    def test_invalid_symbol(self):
        with pytest.raises(ValueError, match="Cannot parse nuclide symbol"):
            _parse_nuclide_symbol("99Tc")

    def test_unknown_element(self):
        with pytest.raises(ValueError, match="Unknown element symbol"):
            _parse_nuclide_symbol("Xx99")


# ---------------------------------------------------------------------------
# validate_json_input
# ---------------------------------------------------------------------------


class TestValidateJsonInput:
    """Test JSON input validation."""

    def test_valid_input(self):
        data = {
            "endpoint_MeV": 0.5,
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
        }
        validate_json_input(data)  # should not raise

    def test_missing_endpoint(self):
        data = {"Z_parent": 43, "Z_daughter": 44, "A_number": 99}
        with pytest.raises(ValueError, match="Missing required field"):
            validate_json_input(data)

    def test_negative_endpoint(self):
        data = {
            "endpoint_MeV": -0.5,
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
        }
        with pytest.raises(ValueError, match="endpoint_MeV must be positive"):
            validate_json_input(data)

    def test_invalid_transition_type(self):
        data = {
            "endpoint_MeV": 0.5,
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
            "transition_type": "invalid",
        }
        with pytest.raises(ValueError, match="Invalid transition_type"):
            validate_json_input(data)


# ---------------------------------------------------------------------------
# load_json_input / json_to_config
# ---------------------------------------------------------------------------


class TestJsonInput:
    """Test JSON file loading and conversion."""

    def _write_json(self, data: dict) -> Path:
        """Write data to a temp JSON file."""
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(data, tmp)
        tmp.close()
        return Path(tmp.name)

    def test_load_and_convert(self):
        data = {
            "parent_symbol": "Tc99",
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
            "endpoint_MeV": 0.294,
            "transition_type": "F2",
            "e_step_MeV": 0.001,
        }
        path = self._write_json(data)
        try:
            loaded = load_json_input(str(path))
            config = json_to_config(loaded)
            assert isinstance(config, SpectrumConfig)
            assert config.Z_parent == 43
            assert config.Z_daughter == 44
            assert config.A_number == 99
            assert abs(config.endpoint_MeV - 0.294) < 1e-6
            assert config.transition_type == "F2"
            assert config.e_step_MeV == 0.001
        finally:
            path.unlink()

    def test_defaults_filled(self):
        data = {
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
            "endpoint_MeV": 0.294,
        }
        path = self._write_json(data)
        try:
            loaded = load_json_input(str(path))
            # Defaults should be filled in
            assert loaded["detector_model"] == "gaussian"
            assert loaded["detector_sigma_a_keV"] == 1.0
        finally:
            path.unlink()

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_json_input("/nonexistent/file.json")


# ---------------------------------------------------------------------------
# decay_info_to_config
# ---------------------------------------------------------------------------


class TestDecayInfoToConfig:
    """Test conversion from DecayInfo to SpectrumConfig."""

    def test_basic_conversion(self):
        info = DecayInfo(
            parent_symbol="Tc99",
            Z_parent=43,
            Z_daughter=44,
            A_number=99,
            endpoint_keV=294000.0,
            endpoint_MeV=0.294,
            transition_type="A",
            forbiddenness_code="0A",
            parent_J=4.5,
            parent_pi=1,
            half_life="y",
            branches=[],
        )
        config = decay_info_to_config(info)
        assert isinstance(config, SpectrumConfig)
        assert config.Z_parent == 43
        assert config.Z_daughter == 44
        assert config.A_number == 99
        assert abs(config.endpoint_MeV - 0.294) < 1e-6
        assert config.transition_type == "A"

    def test_with_detector_response(self):
        info = DecayInfo(
            parent_symbol="Co60",
            Z_parent=27,
            Z_daughter=28,
            A_number=60,
            endpoint_keV=317860.0,
            endpoint_MeV=0.31786,
            transition_type="A",
            forbiddenness_code="1A",
            parent_J=0.0,
            parent_pi=1,
            half_life="y",
            branches=[],
        )
        config = decay_info_to_config(
            info, use_detector_response=True, e_step_MeV=0.005
        )
        assert config.use_detector_response is True
        assert config.e_step_MeV == 0.005


# ---------------------------------------------------------------------------
# get_decay_info_from_paceENSDF
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _has_paceENSDF(), reason="paceENSDF not installed")
class TestPaceENSDFIntegration:
    """Test paceENSDF data retrieval."""

    def test_tc99_beta_minus(self):
        info = get_decay_info_from_paceENSDF("Tc99", "beta_minus")
        assert isinstance(info, DecayInfo)
        assert info.parent_symbol == "Tc99"
        assert info.Z_parent == 43
        assert info.Z_daughter == 44
        assert info.A_number == 99
        assert info.endpoint_MeV > 0
        assert len(info.branches) > 0

    def test_co60_beta_minus(self):
        info = get_decay_info_from_paceENSDF("Co60", "beta_minus")
        assert isinstance(info, DecayInfo)
        assert info.parent_symbol == "Co60"
        assert info.Z_parent == 27
        assert info.Z_daughter == 28
        assert info.A_number == 60
        assert info.endpoint_MeV > 0

    def test_branch_info(self):
        info = get_decay_info_from_paceENSDF("Tc99", "beta_minus")
        for branch in info.branches:
            assert isinstance(branch, BranchInfo)
            assert branch.level_index >= 0
            assert branch.level_energy_keV >= 0
            assert 0 < branch.intensity <= 100
            if branch.log_ft is not None:
                assert branch.log_ft > 0

    def test_invalid_nuclide(self):
        with pytest.raises(ValueError, match="Unknown element symbol"):
            get_decay_info_from_paceENSDF("Xx99", "beta_minus")

    def test_no_decay_data(self):
        # Pb208 is stable, no beta-minus decay data
        with pytest.raises(ValueError, match="No BM decay data found"):
            get_decay_info_from_paceENSDF("Pb208", "beta_minus")


# ---------------------------------------------------------------------------
# create_config_from_source
# ---------------------------------------------------------------------------


class TestCreateConfigFromSource:
    """Test the unified config creation interface."""

    def _write_json(self, data: dict) -> Path:
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(data, tmp)
        tmp.close()
        return Path(tmp.name)

    def test_json_source(self):
        data = {
            "Z_parent": 43,
            "Z_daughter": 44,
            "A_number": 99,
            "endpoint_MeV": 0.294,
        }
        path = self._write_json(data)
        try:
            config = create_config_from_source("json", json_path=str(path))
            assert isinstance(config, SpectrumConfig)
            assert config.Z_parent == 43
        finally:
            path.unlink()

    def test_unknown_source(self):
        with pytest.raises(ValueError, match="Unknown source"):
            create_config_from_source("unknown", nuclide="Tc99")

    def test_paceENSDF_requires_nuclide(self):
        with pytest.raises(ValueError, match="nuclide parameter required"):
            create_config_from_source("paceENSDF")

    def test_json_requires_path(self):
        with pytest.raises(ValueError, match="json_path parameter required"):
            create_config_from_source("json")
