"""Tests for the PSGE-II validation module.

Covers:
- OracleV11 analytical reference values (psge.validation.oracles)
- OracleV12 (raises NotImplementedError — v1.2 in development)
- ValidationSuite (raises NotImplementedError — in development)
"""

import numpy as np
import pytest

from psge.validation.oracles import OracleV11, OracleV12
from psge.validation.suite import ValidationSuite


# ===========================================================================
# OracleV11
# ===========================================================================

class TestOracleV11:
    """Tests for the v1.1 analytical oracle values."""

    def test_flat_square_deficit_is_zero(self):
        """Flat geometry oracle: deficit must be exactly 0."""
        deficit = OracleV11.flat_square_deficit()
        assert deficit == 0.0

    def test_flat_square_deficit_type(self):
        """flat_square_deficit returns a numeric value."""
        assert isinstance(OracleV11.flat_square_deficit(), (int, float))

    def test_cone_deficit_default_value(self):
        """Cone deficit with default parameter = π/2."""
        deficit = OracleV11.cone_deficit()
        assert np.isclose(deficit, np.pi / 2.0)

    def test_cone_deficit_custom_value(self):
        """cone_deficit returns the supplied deficit_value unchanged."""
        custom = np.pi / 3.0
        assert np.isclose(OracleV11.cone_deficit(custom), custom)

    def test_cone_deficit_zero(self):
        """cone_deficit of 0 returns 0 (flat cone)."""
        assert np.isclose(OracleV11.cone_deficit(0.0), 0.0)

    def test_scaling_invariance_returns_dict(self):
        """scaling_invariance returns a dict with expected keys."""
        result = OracleV11.scaling_invariance(2.0)
        assert isinstance(result, dict)
        assert "deficit_invariant" in result
        assert "volume_scales_by" in result

    def test_scaling_invariance_deficit_invariant(self):
        """Deficits are invariant under scaling."""
        result = OracleV11.scaling_invariance(5.0)
        assert result["deficit_invariant"] is True

    def test_scaling_invariance_volume_factor(self):
        """Volumes scale as scale_factor^2 according to the oracle."""
        k = 3.0
        result = OracleV11.scaling_invariance(k)
        assert np.isclose(result["volume_scales_by"], k**2)


# ===========================================================================
# OracleV12 — NotImplementedError (v1.2 in development)
# ===========================================================================

class TestOracleV12:
    """Tests for the v1.2 oracle stubs (currently NotImplementedError)."""

    def test_intrinsic_dihedral_flat_not_implemented(self):
        """intrinsic_dihedral_flat_simplex raises NotImplementedError."""
        with pytest.raises(NotImplementedError):
            OracleV12.intrinsic_dihedral_flat_simplex()

    def test_curved_cone_oracle_not_implemented(self):
        """curved_cone_oracle raises NotImplementedError."""
        with pytest.raises(NotImplementedError):
            OracleV12.curved_cone_oracle(np.pi / 4.0)


# ===========================================================================
# ValidationSuite — NotImplementedError (in development)
# ===========================================================================

class TestValidationSuite:
    """Tests for psge.validation.suite.ValidationSuite."""

    def test_init_v11(self):
        """ValidationSuite can be instantiated with version '1.1'."""
        suite = ValidationSuite(version="1.1")
        assert suite.version == "1.1"

    def test_init_v12(self):
        """ValidationSuite can be instantiated with version '1.2'."""
        suite = ValidationSuite(version="1.2")
        assert suite.version == "1.2"

    def test_run_all_v11_not_implemented(self):
        """ValidationSuite.run_all('1.1') raises NotImplementedError."""
        suite = ValidationSuite(version="1.1")
        with pytest.raises(NotImplementedError):
            suite.run_all()

    def test_run_all_v12_not_implemented(self):
        """ValidationSuite.run_all('1.2') raises NotImplementedError."""
        suite = ValidationSuite(version="1.2")
        with pytest.raises(NotImplementedError):
            suite.run_all()

    def test_run_all_unknown_version_raises_value_error(self):
        """ValidationSuite.run_all raises ValueError for unknown version."""
        suite = ValidationSuite(version="9.9")
        with pytest.raises(ValueError):
            suite.run_all()

    def test_generate_report_not_implemented(self, tmp_path):
        """generate_report raises NotImplementedError (in development)."""
        suite = ValidationSuite(version="1.1")
        with pytest.raises(NotImplementedError):
            suite.generate_report(str(tmp_path / "report.pdf"))
