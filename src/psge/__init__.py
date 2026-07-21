"""PSGE-II: Intrinsic Euclidean Regge Geometry Engine."""

__version__ = "1.2.0-dev"
__author__ = "Jules Demangeot"

from . import core, curvature, mesh, validation

__all__ = ["core", "mesh", "curvature", "validation"]
