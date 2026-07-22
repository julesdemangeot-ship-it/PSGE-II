"""Shared pytest fixtures for PSGE-II test suite."""

import math

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Regular tetrahedron (edge length = 1)
# ---------------------------------------------------------------------------

@pytest.fixture
def regular_tetrahedron_points() -> np.ndarray:
    """Vertices of a regular tetrahedron with unit edge length.

    Returns
    -------
    np.ndarray
        Array of shape (4, 3) containing the four vertex coordinates.
    """
    s3 = math.sqrt(3)
    s6 = math.sqrt(6)
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, s3 / 2.0, 0.0],
            [0.5, s3 / 6.0, s6 / 3.0],
        ]
    )


@pytest.fixture
def regular_tetrahedron_distances() -> np.ndarray:
    """Pairwise distance matrix for a regular tetrahedron (all edges = 1).

    Returns
    -------
    np.ndarray
        Symmetric distance matrix of shape (4, 4).
    """
    d = np.ones((4, 4))
    np.fill_diagonal(d, 0.0)
    return d


@pytest.fixture
def regular_tetrahedron_dihedral_angle() -> float:
    """Exact dihedral angle of a regular tetrahedron.

    Returns
    -------
    float
        arccos(1/3) in radians (~1.2310 rad / ~70.53°).
    """
    return math.acos(1.0 / 3.0)


# ---------------------------------------------------------------------------
# Degenerate simplex (all points collinear → zero volume)
# ---------------------------------------------------------------------------

@pytest.fixture
def degenerate_simplex_points() -> np.ndarray:
    """Four collinear points forming a degenerate (flat) simplex.

    Returns
    -------
    np.ndarray
        Array of shape (4, 3).
    """
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ]
    )


# ---------------------------------------------------------------------------
# Flat triangular pair (two co-planar triangles sharing an edge)
# ---------------------------------------------------------------------------

@pytest.fixture
def flat_dihedral_points() -> tuple:
    """Four points forming two co-planar triangles with a shared edge.

    The shared edge is p1–p2 (along the x-axis).  The two far vertices
    p3 and p4 lie on opposite sides of this edge in the xy-plane, so the
    dihedral angle is π (flat geometry).

    Returns
    -------
    tuple of np.ndarray
        (p1, p2, p3, p4)
    """
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([1.0, 0.0, 0.0])
    p3 = np.array([0.0, 1.0, 0.0])
    p4 = np.array([0.0, -1.0, 0.0])
    return p1, p2, p3, p4


# ---------------------------------------------------------------------------
# Right-angle dihedral (two orthogonal faces)
# ---------------------------------------------------------------------------

@pytest.fixture
def right_angle_dihedral_points() -> tuple:
    """Four points forming a π/2 dihedral angle.

    Shared edge: p1–p2 (x-axis).  p3 is in the xy-plane, p4 in the
    xz-plane → dihedral angle = π/2.

    Returns
    -------
    tuple of np.ndarray
        (p1, p2, p3, p4)
    """
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([1.0, 0.0, 0.0])
    p3 = np.array([0.0, 1.0, 0.0])
    p4 = np.array([0.0, 0.0, 1.0])
    return p1, p2, p3, p4
