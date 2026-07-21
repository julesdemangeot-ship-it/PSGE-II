"""Core geometric and mathematical kernel.

Public API:
- EdgeLengths: Intrinsic metric (edge-length function)
- intrinsic_dihedral: Dihedral angle from edge lengths
- cayley_menger_volume: Volume via Cayley-Menger determinant
- gram_matrix, local_coords, dihedral_from_coords: Low-level utilities
"""

from psge.core.intrinsic_utils import (
    EdgeLengths,
    gram_matrix,
    local_coords,
    dihedral_from_coords,
    cayley_menger_volume,
)
from psge.core.geometry_int import intrinsic_dihedral

__all__ = [
    "EdgeLengths",
    "intrinsic_dihedral",
    "cayley_menger_volume",
    "gram_matrix",
    "local_coords",
    "dihedral_from_coords",
]
