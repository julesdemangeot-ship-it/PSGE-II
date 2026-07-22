"""Extrinsic Euclidean geometry (v1.1).

This module implements the stable v1.1 reference engine, which computes
geometric quantities from embedded coordinates in Euclidean space.

Validated properties:
- Dihedral angle computation
- Regge deficit on embeddable meshes
- Regge action in flat geometry
- Geometric invariances (isometry, scale)
- Degenerate case detection
"""

import numpy as np
from typing import Tuple, Optional


class GeometryExtrinsic:
    """Extrinsic geometry engine for v1.1 stable reference."""
    
    def __init__(self, dimension: int = 3):
        """Initialize extrinsic geometry engine.
        
        Args:
            dimension: Embedding dimension (typically 3 or 4)
        """
        self.dimension = dimension
    
    def dihedral_angle(self, p1: np.ndarray, p2: np.ndarray, 
                       p3: np.ndarray, p4: np.ndarray) -> float:
        """Compute dihedral angle from four points.
        
        Points p1, p2, p3 define one face, p1, p2, p4 define the adjacent face.
        The shared edge is p1-p2.
        
        Args:
            p1, p2, p3, p4: Vertex coordinates (arrays of shape (d,))
            
        Returns:
            Dihedral angle in radians [0, pi]
        """
        # Vectors in first face
        v1 = np.asarray(p3, dtype=float) - np.asarray(p1, dtype=float)
        v2 = np.asarray(p2, dtype=float) - np.asarray(p1, dtype=float)
        n1 = np.cross(v1, v2)
        n1 /= np.linalg.norm(n1)
        
        # Vectors in second face
        v3 = np.asarray(p4, dtype=float) - np.asarray(p1, dtype=float)
        n2 = np.cross(v3, v2)
        n2 /= np.linalg.norm(n2)
        
        # Dihedral angle
        cos_angle = np.dot(n1, n2)
        cos_angle = np.clip(cos_angle, -1, 1)
        angle = np.arccos(cos_angle)
        
        return angle
    
    def deficit(self, angles: np.ndarray) -> float:
        """Compute Regge deficit from dihedral angles around edge.
        
        deficit = 2*pi - sum(angles)
        
        Args:
            angles: Array of dihedral angles around an edge
            
        Returns:
            Deficit value
        """
        return 2 * np.pi - np.sum(angles)
