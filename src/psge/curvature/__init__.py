"""Curvature analysis and Regge action computations.

Public API:
- intrinsic_deficit: Regge deficit at a hinge
- regge_action: Total Regge action
- is_interior_hinge: Topology detection
- hinges_and_cells: Extract hinges and incident cells
"""

from psge.curvature.deficit_int import (
    intrinsic_deficit,
    regge_action,
    is_interior_hinge,
    hinges_and_cells,
)

__all__ = [
    "intrinsic_deficit",
    "regge_action",
    "is_interior_hinge",
    "hinges_and_cells",
]
