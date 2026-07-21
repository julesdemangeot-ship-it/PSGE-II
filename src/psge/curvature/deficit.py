"""Regge deficit computation (extrinsic and intrinsic)."""

import numpy as np
from typing import Dict, Optional


def deficit_extrinsic(dihedral_angles: np.ndarray) -> float:
    """Compute Regge deficit from dihedral angles (v1.1).
    
    deficit = 2*pi - sum(angles)
    
    Args:
        dihedral_angles: Array of angles around an edge
        
    Returns:
        Deficit value
    """
    return 2 * np.pi - np.sum(dihedral_angles)


def deficit_intrinsic(edge_data: Dict) -> Optional[float]:
    """Compute Regge deficit from intrinsic edge lengths (v1.2).
    
    Args:
        edge_data: Dictionary with edge length information
        
    Returns:
        Deficit value or None if degenerate
    """
    # TODO: Implement intrinsic deficit
    raise NotImplementedError("Intrinsic deficit computation in v1.2 development")
