"""Analytical oracles for v1.2 validation.

Provides reference geometries and analytical values for testing the
intrinsic engine:
- Regular tetrahedra (all edges unit length)
- Regular cones (m regular tets around a central edge)
- Flat fans (embeddable reference meshes)

These are the ground truth against which numerical implementations are checked.
"""

from __future__ import annotations

import numpy as np
from itertools import combinations
from typing import Dict, Sequence, Tuple

from psge.core.intrinsic_utils import Edge, Tet, Vertex, EdgeLengths

__all__ = [
    "cone_deficit",
    "regular_cone",
    "regular_tetrahedron",
    "flat_fan",
]


def cone_deficit(m: int) -> float:
    """Analytical deficit of m regular tetrahedra glued around a common edge.

    For a closed ring of m regular unit tetrahedra sharing a central edge (0,1),
    the dihedral angle at each is arccos(1/3). The deficit is:

        deficit = 2π - m * arccos(1/3)

    For m = 6: deficit ≈ 0 (flat, embeddable in R^3)
    For m ≠ 6: deficit ≠ 0 (curved, non-embeddable)

    Args:
        m: Number of tetrahedra in the ring (m >= 3)

    Returns:
        Analytical deficit value
    """
    if m < 3:
        raise ValueError("cone needs at least 3 tetrahedra")
    return float(2.0 * np.pi - m * np.arccos(1.0 / 3.0))


def regular_cone(m: int) -> Tuple[EdgeLengths, Sequence[Tet], Edge]:
    """Closed ring of m regular unit tetrahedra sharing a central edge.

    All edges have length 1, so every cell is a regular tetrahedron.
    For m ≠ 6, this configuration has no isometric embedding in flat R^3 --
    it is intrinsically curved. The central edge (0, 1) is topologically interior.

    Args:
        m: Number of tetrahedra (m >= 3)

    Returns:
        (metric, tetrahedra, central_hinge) where:
        - metric: EdgeLengths with all edges = 1.0
        - tetrahedra: list of m Tet tuples
        - central_hinge: (0, 1)
    """
    if m < 3:
        raise ValueError("cone needs at least 3 tetrahedra")

    A, B = 0, 1
    ring = [2 + k for k in range(m)]
    tets = [(A, B, ring[k], ring[(k + 1) % m]) for k in range(m)]

    # Extract all edges
    edges = set()
    for t in tets:
        for e in combinations(t, 2):
            edges.add(tuple(sorted(e)))

    # All edges have unit length
    metric = EdgeLengths({e: 1.0 for e in edges})

    return metric, tets, (A, B)


def regular_tetrahedron() -> Tuple[Dict[Vertex, np.ndarray], Tet]:
    """Regular unit tetrahedron in global coordinates.

    All edges have length 1. Returns the vertex coordinates and the
    tetrahedron tuple.

    Returns:
        (coords, tet) where:
        - coords: dict mapping vertex indices to 3D coordinates
        - tet: (0, 1, 2, 3)
    """
    # Standard embedding of regular tetrahedron
    a = 1.0 / np.sqrt(2.0)
    coords = {
        0: np.array([0.0, 0.0, 0.0]),
        1: np.array([1.0, 0.0, 0.0]),
        2: np.array([0.5, a, 0.0]),
        3: np.array([0.5, a / 3.0, np.sqrt(2.0 / 3.0)]),
    }
    tet = (0, 1, 2, 3)
    return coords, tet


def flat_fan() -> Tuple[Dict[Vertex, np.ndarray], Sequence[Tet]]:
    """Fan of tetrahedra arranged around a central edge (0, 1) in flat R^3.

    This is an embeddable configuration with zero deficit at the central edge.
    Used for regression testing: intrinsic engine should reproduce zero deficit
    just like the extrinsic engine.

    Returns:
        (coords, tets) where:
        - coords: vertex coordinates in 3D
        - tets: sequence of tetrahedra
    """
    # Central edge along z-axis
    coords = {
        0: np.array([0.0, 0.0, 0.0]),
        1: np.array([0.0, 0.0, 1.0]),
    }

    # Six tetrahedra arranged around the central edge (0, 1)
    # Each has two outer vertices in a ring at z=0.5
    angle_step = 2.0 * np.pi / 6.0
    R = 0.5  # radius of outer ring
    for k in range(6):
        angle = k * angle_step
        idx = 2 + k
        coords[idx] = np.array([R * np.cos(angle), R * np.sin(angle), 0.5])

    ring = [2 + k for k in range(6)]
    tets = [(0, 1, ring[k], ring[(k + 1) % 6]) for k in range(6)]

    return coords, tets
