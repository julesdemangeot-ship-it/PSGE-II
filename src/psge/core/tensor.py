"""Tensor operations for metric computations.

Supports:
- Gram matrices for embedded simplices
- Cayley-Menger determinants for intrinsic geometry
- Metric signatures (Euclidean, Lorentzian)
"""

import numpy as np
from typing import Tuple


def gram_matrix(points: np.ndarray) -> np.ndarray:
    """Compute Gram matrix from a set of points.
    
    Args:
        points: Array of shape (n, d) where n is number of points, d is dimension
        
    Returns:
        Gram matrix of shape (n, n) with G[i,j] = <p_i - p_0, p_j - p_0>
    """
    # Center at first point
    centered = points - points[0]
    return np.dot(centered, centered.T)


def cayley_menger_determinant(distances: np.ndarray) -> float:
    """Compute Cayley-Menger determinant from pairwise distances.
    
    Args:
        distances: Pairwise distance matrix of shape (n, n)
        
    Returns:
        Cayley-Menger determinant value
    """
    n = distances.shape[0]
    # Construct Cayley-Menger matrix
    cm = np.zeros((n + 1, n + 1))
    cm[0, 1:] = 1
    cm[1:, 0] = 1
    cm[1:, 1:] = distances ** 2
    
    return np.linalg.det(cm)


def metric_signature(gram: np.ndarray) -> Tuple[int, int, int]:
    """Compute signature (p, q, z) of metric.
    
    Args:
        gram: Gram matrix or metric tensor
        
    Returns:
        (p, q, z) where p = positive eigenvalues, q = negative, z = zero
    """
    eigenvalues = np.linalg.eigvalsh(gram)
    tolerance = 1e-10
    
    positive = np.sum(eigenvalues > tolerance)
    negative = np.sum(eigenvalues < -tolerance)
    zero = np.sum(np.abs(eigenvalues) <= tolerance)
    
    return positive, negative, zero
