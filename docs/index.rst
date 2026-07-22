PSGE-II Documentation
=====================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   tutorial

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api_reference
   theory

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing

Welcome to PSGE-II
------------------

**PSGE-II** (Intrinsic Euclidean Regge Geometry Engine) is a Python library for
discrete spacetime geometry using Regge calculus.  It provides:

* **Gram matrices** and **Cayley-Menger determinants** for intrinsic computations.
* **Simplex volumes** from both embedded coordinates (extrinsic, v1.1) and edge
  lengths alone (intrinsic, v1.2).
* **Dihedral angles** and **Regge deficits** for curved-geometry analysis.
* **Regge action** — a discrete analogue of the Einstein-Hilbert action.
* **Validation oracles** with known analytical reference values.

Quick Example
-------------

.. code-block:: python

   import numpy as np
   from psge.core.tensor import gram_matrix, cayley_menger_determinant
   from psge.core.volume import simplex_volume_extrinsic
   from psge.core.geometry_ext import GeometryExtrinsic

   # Regular tetrahedron with unit edge length
   import math
   s3, s6 = math.sqrt(3), math.sqrt(6)
   vertices = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.5, s3 / 2, 0.0],
       [0.5, s3 / 6, s6 / 3],
   ])

   # Extrinsic volume
   vol = simplex_volume_extrinsic(vertices)
   print(f"Volume = {vol:.6f}")   # sqrt(2)/12 ≈ 0.117851

   # Dihedral angle
   geo = GeometryExtrinsic()
   angle = geo.dihedral_angle(*vertices[:4])
   print(f"Dihedral = {math.degrees(angle):.2f}°")  # ≈ 70.53°

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
