"""Tests for the PSGE-II core geometric kernel.

Covers:
- Gram matrices (psge.core.tensor)
- Cayley-Menger determinants (psge.core.tensor)
- Metric signatures (psge.core.tensor)
- Simplex volumes, extrinsic and intrinsic (psge.core.volume)
- Dihedral angles and deficits (psge.core.geometry_ext)
"""

import math

import numpy as np
import pytest

from psge.core.tensor import cayley_menger_determinant, gram_matrix, metric_signature
from psge.core.volume import simplex_volume_extrinsic, simplex_volume_intrinsic
from psge.core.geometry_ext import GeometryExtrinsic


# ===========================================================================
# Gram matrix
# ===========================================================================

class TestGramMatrix:
    """Tests for gram_matrix."""

    def test_gram_matrix_shape(self, regular_tetrahedron_points):
        """Gram matrix should be square with size equal to number of points."""
        pts = regular_tetrahedron_points
        G = gram_matrix(pts)
        n = pts.shape[0]
        assert G.shape == (n, n)

    def test_gram_matrix_symmetry(self, regular_tetrahedron_points):
        """Gram matrix must be symmetric."""
        G = gram_matrix(regular_tetrahedron_points)
        assert np.allclose(G, G.T)

    def test_gram_matrix_first_row_zero(self, regular_tetrahedron_points):
        """First row and column of Gram matrix are zero (centered at p0)."""
        G = gram_matrix(regular_tetrahedron_points)
        assert np.allclose(G[0, :], 0.0)
        assert np.allclose(G[:, 0], 0.0)

    def test_gram_matrix_regular_tetrahedron(self, regular_tetrahedron_points):
        """Gram matrix of a regular tetrahedron has the expected structure.

        For a unit regular tetrahedron centred at p0:
        - Diagonal entries (except [0,0]) should equal 1  (|e_i|^2 = 1)
        - Off-diagonal entries (i,j both > 0) should equal 0.5 (cos(60°))
        """
        G = gram_matrix(regular_tetrahedron_points)
        # Diagonal entries (skip index 0 which is always 0)
        for i in range(1, 4):
            assert np.isclose(G[i, i], 1.0, atol=1e-10), f"G[{i},{i}] = {G[i,i]}"
        # Off-diagonal
        for i in range(1, 4):
            for j in range(1, 4):
                if i != j:
                    assert np.isclose(G[i, j], 0.5, atol=1e-10), f"G[{i},{j}] = {G[i,j]}"

    def test_gram_matrix_single_point(self):
        """Gram matrix of a single point is the 1×1 zero matrix."""
        G = gram_matrix(np.array([[3.0, 4.0, 5.0]]))
        assert G.shape == (1, 1)
        assert G[0, 0] == 0.0


# ===========================================================================
# Cayley-Menger determinant
# ===========================================================================

class TestCayleyMengerDeterminant:
    """Tests for cayley_menger_determinant."""

    def test_degenerate_distance_matrix(self):
        """Three collinear points have a zero (degenerate) CM determinant."""
        # Three equally-spaced collinear points
        d = np.array([
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [2.0, 1.0, 0.0],
        ])
        cm = cayley_menger_determinant(d)
        assert np.isclose(cm, 0.0, atol=1e-10)

    def test_equilateral_triangle(self):
        """CM determinant for an equilateral triangle with side 1.

        For a 2-simplex with all sides = 1:
        (n! * V)^2 = (-1)^(n+1) / 2^n * det(CM)
        For n=2, V = sqrt(3)/4:
        (2! * sqrt(3)/4)^2 = (-1)^3 / 4 * det(CM)
        (sqrt(3)/2)^2 = -det(CM)/4
        det(CM) = -3
        """
        d = np.array([
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ])
        cm = cayley_menger_determinant(d)
        assert np.isclose(cm, -3.0, atol=1e-10)

    def test_returns_float(self, regular_tetrahedron_distances):
        """Return value should be a float-compatible scalar."""
        cm = cayley_menger_determinant(regular_tetrahedron_distances)
        assert np.isfinite(cm)


# ===========================================================================
# Metric signature
# ===========================================================================

class TestMetricSignature:
    """Tests for metric_signature."""

    def test_identity_matrix_is_euclidean(self):
        """3×3 identity: signature (3, 0, 0)."""
        p, q, z = metric_signature(np.eye(3))
        assert (p, q, z) == (3, 0, 0)

    def test_zero_matrix_signature(self):
        """3×3 zero matrix: signature (0, 0, 3)."""
        p, q, z = metric_signature(np.zeros((3, 3)))
        assert (p, q, z) == (0, 0, 3)

    def test_regular_tetrahedron_gram_signature(self, regular_tetrahedron_points):
        """Gram matrix of a regular tetrahedron is PSD with one zero eigenvalue."""
        G = gram_matrix(regular_tetrahedron_points)
        p, q, z = metric_signature(G)
        # Centered at p0, so one zero; others positive
        assert q == 0
        assert z >= 1


# ===========================================================================
# Simplex volume — extrinsic
# ===========================================================================

class TestSimplexVolumeExtrinsic:
    """Tests for simplex_volume_extrinsic."""

    def test_unit_cube_tetrahedron_volume(self):
        """Tetrahedron inscribed in the unit cube: V = 1/6."""
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        vol = simplex_volume_extrinsic(points)
        assert np.isclose(vol, 1.0 / 6.0, atol=1e-12)

    def test_regular_tetrahedron_volume(self, regular_tetrahedron_points):
        """Regular tetrahedron with edge=1 has V = sqrt(2)/12."""
        vol = simplex_volume_extrinsic(regular_tetrahedron_points)
        expected = math.sqrt(2) / 12.0
        assert np.isclose(vol, expected, atol=1e-12)

    def test_degenerate_simplex_zero_volume(self, degenerate_simplex_points):
        """Collinear points produce zero volume."""
        vol = simplex_volume_extrinsic(degenerate_simplex_points)
        assert np.isclose(vol, 0.0, atol=1e-12)

    def test_scaling_behaviour(self, regular_tetrahedron_points):
        """Scaling all coordinates by k scales volume by k^3."""
        k = 3.0
        vol_original = simplex_volume_extrinsic(regular_tetrahedron_points)
        vol_scaled = simplex_volume_extrinsic(regular_tetrahedron_points * k)
        assert np.isclose(vol_scaled, vol_original * k**3, rtol=1e-10)

    def test_translation_invariance(self, regular_tetrahedron_points):
        """Translating all vertices does not change the volume."""
        translation = np.array([10.0, -5.0, 3.0])
        vol_original = simplex_volume_extrinsic(regular_tetrahedron_points)
        vol_translated = simplex_volume_extrinsic(regular_tetrahedron_points + translation)
        assert np.isclose(vol_original, vol_translated, atol=1e-12)


# ===========================================================================
# Simplex volume — intrinsic
# ===========================================================================

class TestSimplexVolumeIntrinsic:
    """Tests for simplex_volume_intrinsic."""

    def test_equilateral_triangle_area(self):
        """Equilateral triangle (side=1): area = sqrt(3)/4."""
        d = np.array([
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ])
        area = simplex_volume_intrinsic(d)
        assert area is not None
        assert np.isclose(area, math.sqrt(3) / 4.0, atol=1e-12)

    def test_regular_tetrahedron_volume(self, regular_tetrahedron_distances):
        """Regular tetrahedron (edge=1): V = sqrt(2)/12."""
        vol = simplex_volume_intrinsic(regular_tetrahedron_distances)
        assert vol is not None
        expected = math.sqrt(2) / 12.0
        assert np.isclose(vol, expected, atol=1e-10)

    def test_degenerate_returns_none(self):
        """Collinear distance matrix returns None (degenerate simplex)."""
        d = np.array([
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.0],
            [2.0, 1.0, 0.0],
        ])
        result = simplex_volume_intrinsic(d)
        assert result is None

    def test_consistency_with_extrinsic(self, regular_tetrahedron_points, regular_tetrahedron_distances):
        """Intrinsic and extrinsic volumes agree for a regular tetrahedron."""
        vol_ext = simplex_volume_extrinsic(regular_tetrahedron_points)
        vol_int = simplex_volume_intrinsic(regular_tetrahedron_distances)
        assert vol_int is not None
        assert np.isclose(vol_ext, vol_int, atol=1e-10)


# ===========================================================================
# GeometryExtrinsic — dihedral angles and deficits
# ===========================================================================

class TestGeometryExtrinsic:
    """Tests for psge.core.geometry_ext.GeometryExtrinsic."""

    @pytest.fixture
    def geo(self):
        return GeometryExtrinsic(dimension=3)

    # -- dihedral_angle -------------------------------------------------------

    def test_dihedral_flat_is_pi(self, geo, flat_dihedral_points):
        """Co-planar faces → dihedral angle = π."""
        p1, p2, p3, p4 = flat_dihedral_points
        angle = geo.dihedral_angle(p1, p2, p3, p4)
        assert np.isclose(angle, np.pi, atol=1e-10)

    def test_dihedral_right_angle(self, geo, right_angle_dihedral_points):
        """Orthogonal faces → dihedral angle = π/2."""
        p1, p2, p3, p4 = right_angle_dihedral_points
        angle = geo.dihedral_angle(p1, p2, p3, p4)
        assert np.isclose(angle, np.pi / 2, atol=1e-10)

    def test_dihedral_range(self, geo, regular_tetrahedron_points):
        """Dihedral angle of regular tetrahedron is in [0, π]."""
        pts = regular_tetrahedron_points
        angle = geo.dihedral_angle(pts[0], pts[1], pts[2], pts[3])
        assert 0.0 <= angle <= np.pi

    def test_dihedral_regular_tetrahedron(self, geo, regular_tetrahedron_points,
                                          regular_tetrahedron_dihedral_angle):
        """Dihedral angle of regular tetrahedron = arccos(1/3)."""
        pts = regular_tetrahedron_points
        angle = geo.dihedral_angle(pts[0], pts[1], pts[2], pts[3])
        assert np.isclose(angle, regular_tetrahedron_dihedral_angle, atol=1e-10)

    def test_dihedral_accepts_float_arrays(self, geo):
        """dihedral_angle works for explicit float64 arrays."""
        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([1.0, 0.0, 0.0])
        p3 = np.array([0.5, 1.0, 0.0])
        p4 = np.array([0.5, -1.0, 0.0])
        angle = geo.dihedral_angle(p1, p2, p3, p4)
        assert np.isfinite(angle)

    def test_dihedral_accepts_integer_arrays(self, geo):
        """dihedral_angle handles integer-dtype input without error."""
        p1 = np.array([0, 0, 0])
        p2 = np.array([1, 0, 0])
        p3 = np.array([0, 1, 0])
        p4 = np.array([0, -1, 0])
        angle = geo.dihedral_angle(p1, p2, p3, p4)
        assert np.isclose(angle, np.pi, atol=1e-10)

    # -- deficit --------------------------------------------------------------

    def test_deficit_flat_mesh(self, geo):
        """Six angles of π/3 around an edge → deficit = 0 (flat geometry)."""
        angles = np.full(6, np.pi / 3)
        assert np.isclose(geo.deficit(angles), 0.0, atol=1e-14)

    def test_deficit_positive_curvature(self, geo):
        """Five angles of π/3 → positive deficit (positive curvature)."""
        angles = np.full(5, np.pi / 3)
        deficit = geo.deficit(angles)
        assert deficit > 0.0

    def test_deficit_negative_curvature(self, geo):
        """Seven angles of π/3 → negative deficit (negative curvature)."""
        angles = np.full(7, np.pi / 3)
        deficit = geo.deficit(angles)
        assert deficit < 0.0

    def test_deficit_formula(self, geo):
        """deficit = 2π − Σ(angles)."""
        angles = np.array([0.5, 0.7, 0.9, 1.1])
        expected = 2 * np.pi - angles.sum()
        assert np.isclose(geo.deficit(angles), expected)
