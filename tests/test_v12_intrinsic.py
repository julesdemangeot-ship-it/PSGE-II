"""Tests for v1.2 intrinsic formulation (in development)."""

import math

import numpy as np
import pytest

from psge.core.geometry_int import GeometryIntrinsic


def _distance_matrix(points: np.ndarray) -> np.ndarray:
    """Build a pairwise distance matrix from an array of point coordinates."""
    n = len(points)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.linalg.norm(points[i] - points[j])
    return D


class TestGeometryIntrinsic:
    """Test suite for intrinsic geometry (v1.2)."""

    @pytest.fixture
    def geo(self) -> GeometryIntrinsic:
        return GeometryIntrinsic()

    # -- dihedral_angle_from_distances ----------------------------------------

    def test_regular_tetrahedron(self, geo):
        """Regular tetrahedron (all edges = 1): dihedral angle = arccos(1/3)."""
        D = np.ones((4, 4))
        np.fill_diagonal(D, 0.0)
        angle = geo.dihedral_angle_from_distances(D)
        assert angle is not None
        assert np.isclose(angle, math.acos(1.0 / 3.0), atol=1e-10)

    def test_right_angle_dihedral(self, geo):
        """Configuration with orthogonal faces yields dihedral angle π/2."""
        # p0=[0,0,0], p1=[1,0,0] (shared edge), p2=[0,1,0], p3=[0,0,1]
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        D = _distance_matrix(points)
        angle = geo.dihedral_angle_from_distances(D)
        assert angle is not None
        assert np.isclose(angle, math.pi / 2.0, atol=1e-10)

    def test_flat_dihedral_is_pi(self, geo):
        """Co-planar configuration yields dihedral angle π."""
        # p3 on opposite side of edge in the same plane as p2
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
        D = _distance_matrix(points)
        angle = geo.dihedral_angle_from_distances(D)
        assert angle is not None
        assert np.isclose(angle, math.pi, atol=1e-10)

    def test_result_in_range(self, geo):
        """Dihedral angle is always in [0, π]."""
        D = np.ones((4, 4))
        np.fill_diagonal(D, 0.0)
        angle = geo.dihedral_angle_from_distances(D)
        assert angle is not None
        assert 0.0 <= angle <= math.pi

    def test_returns_float(self, geo):
        """Return type is a Python float (or None)."""
        D = np.ones((4, 4))
        np.fill_diagonal(D, 0.0)
        angle = geo.dihedral_angle_from_distances(D)
        assert isinstance(angle, float)

    def test_wrong_shape_raises_value_error(self, geo):
        """Non-4×4 distance matrix raises ValueError."""
        with pytest.raises(ValueError):
            geo.dihedral_angle_from_distances(np.ones((3, 3)))

    def test_zero_edge_length_returns_none(self, geo):
        """Distance matrix with zero-length shared edge returns None."""
        D = np.ones((4, 4))
        np.fill_diagonal(D, 0.0)
        # Make vertices 0 and 1 coincide
        D[0, 1] = 0.0
        D[1, 0] = 0.0
        result = geo.dihedral_angle_from_distances(D)
        assert result is None

    def test_consistency_with_extrinsic(self, geo):
        """Intrinsic result agrees with extrinsic computation for a known tetrahedron."""
        from psge.core.geometry_ext import GeometryExtrinsic

        sqrt_3 = math.sqrt(3)
        sqrt_6 = math.sqrt(6)
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, sqrt_3 / 2.0, 0.0],
            [0.5, sqrt_3 / 6.0, sqrt_6 / 3.0],
        ])
        D = _distance_matrix(points)

        geo_ext = GeometryExtrinsic(dimension=3)
        angle_ext = geo_ext.dihedral_angle(*points[:4])
        angle_int = geo.dihedral_angle_from_distances(D)

        assert angle_int is not None
        assert np.isclose(angle_int, angle_ext, atol=1e-10)

    # -- placeholder ----------------------------------------------------------

    def test_placeholder(self):
        """Placeholder kept for backward compatibility."""
        assert True
