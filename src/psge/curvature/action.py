"""Regge action computation."""

import numpy as np
from typing import Dict, List


def regge_action(deficits: np.ndarray, volumes: np.ndarray) -> float:
    """Compute total Regge action.
    
    S = sum(deficit_i * volume_i)
    
    Args:
        deficits: Array of Regge deficits at each edge
        volumes: Array of dual volumes at each edge
        
    Returns:
        Total Regge action value
    """
    return np.sum(deficits * volumes)
