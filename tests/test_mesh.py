"""Tests for the PSGE-II mesh module.

Covers:
- Hinge extraction (psge.mesh.traversal) — currently raises NotImplementedError
- Barycentric dual volume (psge.mesh.dual) — currently raises NotImplementedError
"""

import pytest
import numpy as np

from psge.mesh.traversal import extract_hinges
from psge.mesh.dual import barycentric_dual_volume


# ===========================================================================
# Hinge extraction
# ===========================================================================

class TestExtractHinges:
    """Tests for psge.mesh.traversal.extract_hinges."""

    def test_extract_hinges_not_implemented(self):
        """extract_hinges raises NotImplementedError (v1.2 in development)."""
        mesh = {}
        with pytest.raises(NotImplementedError):
            extract_hinges(mesh)

    def test_extract_hinges_returns_not_implemented_for_any_input(self):
        """extract_hinges raises NotImplementedError regardless of input."""
        with pytest.raises(NotImplementedError):
            extract_hinges({"cells": [], "vertices": []})


# ===========================================================================
# Barycentric dual volume
# ===========================================================================

class TestBarycentricrDualVolume:
    """Tests for psge.mesh.dual.barycentric_dual_volume."""

    def test_barycentric_dual_not_implemented(self, regular_tetrahedron_points):
        """barycentric_dual_volume raises NotImplementedError (v1.2 in development)."""
        with pytest.raises(NotImplementedError):
            barycentric_dual_volume(regular_tetrahedron_points, 0)

    def test_barycentric_dual_any_vertex_index(self, regular_tetrahedron_points):
        """NotImplementedError for every vertex index."""
        for i in range(4):
            with pytest.raises(NotImplementedError):
                barycentric_dual_volume(regular_tetrahedron_points, i)
