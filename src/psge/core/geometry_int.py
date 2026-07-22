"""Intrinsic Euclidean geometry (v1.2 - In Development).

This module implements the intrinsic formulation of PSGE-II v1.2,
computing geometric quantities from edge lengths alone.

Featured developments:
- Dihedral angles from Cayley-Menger determinants
- Intrinsic Regge deficit computation
- Barycentric dual measures
- Curved cone oracle with non-zero deficit
"""

import numpy as np
from typing import Optional, Tuple


class GeometryIntrinsic:
    """Intrinsic geometry engine for v1.2 development."""
    
    def __init__(self):
        """Initialize intrinsic geometry engine."""
        pass
    
    def dihedral_angle_from_distances(self, distances: np.ndarray) -> Optional[float]:
        """Compute dihedral angle at the edge (vertex 0, vertex 1) from a 4-vertex
        pairwise distance matrix using Cayley-Menger theory.

        The two triangular faces sharing the edge are (0, 1, 2) and (0, 1, 3).
        An embedding in ℝ³ is reconstructed from the distance matrix via the
        Gram matrix; the dihedral angle is then derived from face normals.

        Args:
            distances: Symmetric 4×4 distance matrix with zero diagonal.
                       distances[i, j] = length of edge between vertex i and j.

        Returns:
            Dihedral angle in radians [0, π], or None if the configuration is
            degenerate (zero-length edge or collapsed face).

        Raises:
            ValueError: If *distances* does not have shape (4, 4).
        """
        D = np.asarray(distances, dtype=float)
        if D.shape != (4, 4):
            raise ValueError("distances must be a 4×4 pairwise distance matrix")

        # Build the 3×3 Gram matrix centred at vertex 0.
        # G[i, j] = (D[0, i+1]² + D[0, j+1]² − D[i+1, j+1]²) / 2
        G = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                G[i, j] = (D[0, i + 1] ** 2 + D[0, j + 1] ** 2 - D[i + 1, j + 1] ** 2) / 2.0

        # Recover embedding coordinates via eigendecomposition G = V Λ Vᵀ.
        eigenvalues, eigenvectors = np.linalg.eigh(G)

        if np.any(eigenvalues < -1e-10):
            return None  # Non-embeddable configuration

        eigenvalues = np.maximum(eigenvalues, 0.0)  # clip numerical noise
        # Build P = V @ diag(sqrt(Λ)) so that P @ Pᵀ = G.
        # Broadcasting: each column j of eigenvectors is multiplied by
        # sqrt(eigenvalues[j]), which is equivalent to
        # eigenvectors @ np.diag(np.sqrt(eigenvalues)).
        # Row i of P is the embedding coordinate of vertex i+1.
        coords = eigenvectors * np.sqrt(eigenvalues)  # shape (3, 3)

        p0 = np.zeros(3)
        p1 = coords[0]  # vertex 1
        p2 = coords[1]  # vertex 2
        p3 = coords[2]  # vertex 3

        # Edge vector
        edge = p1 - p0
        edge_len = np.linalg.norm(edge)
        if edge_len < 1e-15:
            return None  # Degenerate: edge has zero length

        edge_unit = edge / edge_len

        # Project p2 and p3 onto the plane perpendicular to the edge.
        r2 = (p2 - p0) - np.dot(p2 - p0, edge_unit) * edge_unit
        r3 = (p3 - p0) - np.dot(p3 - p0, edge_unit) * edge_unit

        norm_r2 = np.linalg.norm(r2)
        norm_r3 = np.linalg.norm(r3)
        if norm_r2 < 1e-15 or norm_r3 < 1e-15:
            return None  # Degenerate: a face is collinear with the edge

        cos_angle = np.dot(r2, r3) / (norm_r2 * norm_r3)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return float(np.arccos(cos_angle))
    
    def deficit_intrinsic(self, edge_configuration: dict) -> Optional[float]:
        """Compute intrinsic Regge deficit from edge lengths.
        
        Args:
            edge_configuration: Dictionary of edge lengths
            
        Returns:
            Deficit value or None if invalid configuration
        """
        # TODO: Implement intrinsic deficit computation
        raise NotImplementedError("v1.2 intrinsic deficit computation in development")
