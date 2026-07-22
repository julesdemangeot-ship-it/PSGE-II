"""Volume and simplex computations.

Provides functions for computing volumes of simplices using both
embedding coordinates (extrinsic) and intrinsic distance data.
"""

import math

import numpy as np
from typing import Optional


def simplex_volume_extrinsic(points: np.ndarray) -> float:
    """Compute volume of simplex from embedded coordinates.

    Uses the formula V = abs(det(edges)) / n! where *edges* are the vectors
    from the first vertex to each of the remaining vertices.

    Args:
        points: Array of shape (n+1, d) representing n-simplex vertices in d-space

    Returns:
        Volume of the simplex
    """
    edges = points[1:] - points[0]
    volume = np.abs(np.linalg.det(edges)) / math.factorial(len(edges))
    return volume


def simplex_volume_intrinsic(distances: np.ndarray) -> Optional[float]:
    """Compute volume of simplex from pairwise distances (intrinsic).
    
    Uses Cayley-Menger determinant formula:
    (n! * V)^2 = (-1)^(n+1) / 2^n * det(CM)
    
    Args:
        distances: Pairwise distance matrix of shape (n, n)
        
    Returns:
        Volume if positive definite, None if degenerate
    """
    from .tensor import cayley_menger_determinant
    
    n = distances.shape[0] - 1  # dimension of simplex
    cm_det = cayley_menger_determinant(distances)
    
    # Check if valid (CM determinant should have correct sign)
    expected_sign = (-1) ** (n + 1)
    if cm_det * expected_sign <= 0:
        return None  # Degenerate or non-embeddable
    
    volume_squared = (expected_sign * cm_det) / (2 ** n * (math.factorial(n) ** 2))
    
    if volume_squared < 0:
        return None
    
    return np.sqrt(volume_squared)
