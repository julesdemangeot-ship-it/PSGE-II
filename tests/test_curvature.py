"""Tests for the PSGE-II curvature module.

Covers:
- Extrinsic deficit computation (psge.curvature.deficit)
- Intrinsic deficit (raises NotImplementedError — v1.2 in development)
- Regge action (psge.curvature.action)
"""

import numpy as np
import pytest

from psge.curvature.deficit import deficit_extrinsic, deficit_intrinsic
from psge.curvature.action import regge_action


# ===========================================================================
# deficit_extrinsic
# ===========================================================================

class TestDeficitExtrinsic:
    """Tests for psge.curvature.deficit.deficit_extrinsic."""

    def test_flat_geometry_zero_deficit(self):
        """Six angles of π/3 around an edge → deficit = 0 (flat geometry)."""
        angles = np.full(6, np.pi / 3.0)
        assert np.isclose(deficit_extrinsic(angles), 0.0, atol=1e-14)

    def test_positive_curvature(self):
        """Fewer than six π/3 angles → positive deficit (positive curvature)."""
        angles = np.full(5, np.pi / 3.0)
        assert deficit_extrinsic(angles) > 0.0

    def test_negative_curvature(self):
        """More than six π/3 angles → negative deficit (negative curvature)."""
        angles = np.full(7, np.pi / 3.0)
        assert deficit_extrinsic(angles) < 0.0

    def test_formula_2pi_minus_sum(self):
        """deficit = 2π − Σ(angles)."""
        angles = np.array([0.3, 0.5, 0.7, 0.9, 1.1])
        expected = 2.0 * np.pi - angles.sum()
        assert np.isclose(deficit_extrinsic(angles), expected)

    def test_single_angle_two_pi(self):
        """A single angle equal to 2π produces zero deficit."""
        assert np.isclose(deficit_extrinsic(np.array([2.0 * np.pi])), 0.0, atol=1e-14)

    def test_empty_angles_returns_two_pi(self):
        """No angles → deficit = 2π (no triangles meet at the edge)."""
        assert np.isclose(deficit_extrinsic(np.array([])), 2.0 * np.pi)

    def test_cone_deficit_pi_over_two(self):
        """Cone geometry: 4 angles of π*3/8 → deficit = π/2."""
        cone_angle = np.pi * 3.0 / 8.0  # sum = 3π/2, deficit = π/2
        angles = np.full(4, cone_angle)
        expected = 2.0 * np.pi - 4.0 * cone_angle
        assert np.isclose(deficit_extrinsic(angles), expected)


# ===========================================================================
# deficit_intrinsic — NotImplementedError (v1.2 in development)
# ===========================================================================

class TestDeficitIntrinsic:
    """Tests for psge.curvature.deficit.deficit_intrinsic."""

    def test_not_implemented(self):
        """deficit_intrinsic raises NotImplementedError (v1.2 in development)."""
        with pytest.raises(NotImplementedError):
            deficit_intrinsic({})

    def test_not_implemented_for_non_empty_dict(self):
        """NotImplementedError for any edge_data input."""
        with pytest.raises(NotImplementedError):
            deficit_intrinsic({"edge": 1.0})


# ===========================================================================
# regge_action
# ===========================================================================

class TestReggeAction:
    """Tests for psge.curvature.action.regge_action."""

    def test_flat_geometry_zero_action(self):
        """All-zero deficits → Regge action = 0 (flat geometry)."""
        deficits = np.zeros(5)
        volumes = np.ones(5)
        assert np.isclose(regge_action(deficits, volumes), 0.0)

    def test_single_edge(self):
        """S = deficit × volume for a single edge."""
        deficits = np.array([np.pi / 4.0])
        volumes = np.array([2.0])
        assert np.isclose(regge_action(deficits, volumes), np.pi / 2.0)

    def test_multiple_edges(self):
        """S = Σ(deficit_i × volume_i)."""
        deficits = np.array([0.1, 0.2, 0.3])
        volumes = np.array([1.0, 2.0, 3.0])
        expected = 0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0  # 1.4
        assert np.isclose(regge_action(deficits, volumes), expected)

    def test_scaling_by_volume(self):
        """Doubling all dual volumes doubles the action."""
        deficits = np.array([0.5, 0.3, 0.7])
        volumes = np.array([1.0, 2.0, 0.5])
        S1 = regge_action(deficits, volumes)
        S2 = regge_action(deficits, 2.0 * volumes)
        assert np.isclose(S2, 2.0 * S1)

    def test_scaling_by_deficit(self):
        """Doubling all deficits doubles the action."""
        deficits = np.array([0.5, 0.3, 0.7])
        volumes = np.array([1.0, 2.0, 0.5])
        S1 = regge_action(deficits, volumes)
        S2 = regge_action(2.0 * deficits, volumes)
        assert np.isclose(S2, 2.0 * S1)

    def test_returns_float(self):
        """Return value should be a scalar float."""
        result = regge_action(np.array([0.1, 0.2]), np.array([1.0, 1.0]))
        assert np.isscalar(result) or result.ndim == 0
