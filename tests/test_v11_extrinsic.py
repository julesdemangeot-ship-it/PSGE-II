"""Tests for v1.1 extrinsic Euclidean engine."""

import pytest
import numpy as np
from psge.core.geometry_ext import GeometryExtrinsic


class TestGeometryExtrinsic:
    """Test suite for extrinsic geometry."""
    
    @pytest.fixture
    def geometry(self):
        """Create geometry engine instance."""
        return GeometryExtrinsic(dimension=3)
    
    def test_dihedral_angle_flat(self, geometry):
        """Test dihedral angle computation for flat configuration."""
        # Four coplanar points
        p1 = np.array([0, 0, 0])
        p2 = np.array([1, 0, 0])
        p3 = np.array([0, 1, 0])
        p4 = np.array([0, -1, 0])
        
        angle = geometry.dihedral_angle(p1, p2, p3, p4)
        
        # For flat configuration, dihedral should be pi
        assert np.isclose(angle, np.pi, atol=1e-10)
    
    def test_deficit_flat_mesh(self, geometry):
        """Test deficit computation for flat mesh."""
        # Six angles around an edge in flat mesh
        angles = np.array([np.pi/3] * 6)
        
        deficit = geometry.deficit(angles)
        
        # For flat mesh, deficit should be 0
        assert np.isclose(deficit, 0, atol=1e-10)
