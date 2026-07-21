"""Intrinsic Regge deficit and action computation.

Computes the Regge deficit at hinges and the total Regge action from
intrinsic edge-length data alone. The deficit measures curvature: zero
on flat regions, non-zero on curved simplicial complexes.

Unlike the extrinsic engine (v1.1), interior hinges can have non-zero
deficit -- this is the classical Regge construction.
"""

from __future__ import annotations

import numpy as np
from typing import Sequence, Tuple, Set, Dict
from itertools import combinations

from psge.core.intrinsic_utils import Edge, Tet, EdgeLengths
from psge.core.geometry_int import intrinsic_dihedral

Vertex = int


def hinges_and_cells(tets: Sequence[Tet]) -> Dict[Edge, list]:
    """Extract all hinges and their incident cells.

    Returns a dictionary mapping each hinge (edge) to the list of
    tetrahedra that contain it.
    """
    hinge_cells: Dict[Edge, list] = {}
    for tet in tets:
        for i, j in combinations(range(4), 2):
            hinge = tuple(sorted([tet[i], tet[j]]))
            if hinge not in hinge_cells:
                hinge_cells[hinge] = []
            hinge_cells[hinge].append(tet)
    return hinge_cells


def is_interior_hinge(hinge: Edge, cells: Sequence[Tet]) -> bool:
    """A hinge is interior iff every incident face is shared by two cells.

    Boundary hinges have at least one incident face with only one cell.
    """
    face_count: Dict[Tuple[int, int, int], int] = {}
    for c in cells:
        for f in combinations(sorted(c), 3):
            face_count[f] = face_count.get(f, 0) + 1

    a, b = hinge
    for c in cells:
        for f in combinations(sorted(c), 3):
            if a in f and b in f and face_count[tuple(sorted(f))] == 1:
                return False
    return True


def intrinsic_deficit(hinge: Edge, cells: Sequence[Tet], 
                      metric: EdgeLengths) -> float:
    """Regge deficit at `hinge`: 2π - Σ(incident dihedral angles).

    Unlike extrinsic (v1.1), this can be non-zero at an interior hinge
    because each cell reconstructs in its own local frame.
    
    Args:
        hinge: Edge (pair of vertices)
        cells: Tetrahedra incident to this hinge
        metric: Intrinsic edge-length metric
        
    Returns:
        Deficit value (non-zero indicates curvature)
    """
    total_angle = sum(intrinsic_dihedral(c, hinge, metric) for c in cells)
    return float(2.0 * np.pi - total_angle)


def hinge_volume(hinge: Edge, cells: Sequence[Tet], 
                 metric: EdgeLengths) -> float:
    """Dual volume (barycentric) at a hinge (sum of dual volumes in incident cells).

    For simplicity, use the sum of cell volumes divided by (2*D+2) where D=3.
    This is a placeholder; the full barycentric measure requires more
    sophisticated dual geometry.
    """
    from psge.core.intrinsic_utils import cayley_menger_volume
    
    total_volume = 0.0
    for c in cells:
        vol = cayley_menger_volume(c, metric)
        # Dual volume at hinge: approximately cell_volume / (2*dim+2)
        # For 3D: divide by 8
        total_volume += vol / 8.0
    return total_volume


def regge_action(hinges: Sequence[Edge], cells_per_hinge: Dict[Edge, Sequence[Tet]],
                 metric: EdgeLengths) -> float:
    """Total Regge action: Σ (deficit_i × volume_i) over interior hinges.

    The Regge action is the discrete analog of the Einstein-Hilbert action
    in General Relativity. On flat spacetime it is zero; on curved spaces
    it is non-zero.

    Args:
        hinges: List of edges (hinges) to include
        cells_per_hinge: Dictionary mapping each hinge to incident cells
        metric: Intrinsic edge-length metric

    Returns:
        Total action value
    """
    action = 0.0
    for hinge in hinges:
        cells = cells_per_hinge[hinge]
        if is_interior_hinge(hinge, cells):
            deficit = intrinsic_deficit(hinge, cells, metric)
            volume = hinge_volume(hinge, cells, metric)
            action += deficit * volume
    return float(action)
