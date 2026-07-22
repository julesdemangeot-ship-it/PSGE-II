Tutorial
========

This tutorial walks through the main use-cases of PSGE-II step by step.

.. contents::
   :local:
   :depth: 2

Computing Simplex Volumes
-------------------------

PSGE-II provides two ways to compute volumes: from embedded coordinates
(*extrinsic*) and from edge lengths alone (*intrinsic*).

Extrinsic Volume
~~~~~~~~~~~~~~~~

The extrinsic approach uses explicit vertex coordinates:

.. code-block:: python

   import numpy as np
   from psge.core.volume import simplex_volume_extrinsic

   # Tetrahedron inscribed in the unit cube
   points = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ])

   vol = simplex_volume_extrinsic(points)
   print(f"Volume = {vol:.6f}")  # 1/6 ≈ 0.166667

Intrinsic Volume
~~~~~~~~~~~~~~~~

The intrinsic approach uses only pairwise distances (Cayley-Menger):

.. code-block:: python

   import numpy as np
   from psge.core.volume import simplex_volume_intrinsic

   # Regular tetrahedron: all edges = 1
   d = np.ones((4, 4))
   np.fill_diagonal(d, 0.0)

   vol = simplex_volume_intrinsic(d)
   print(f"Volume = {vol:.6f}")  # sqrt(2)/12 ≈ 0.117851

   # Degenerate simplex (collinear points) returns None
   d_degen = np.array([
       [0., 1., 2.],
       [1., 0., 1.],
       [2., 1., 0.],
   ])
   assert simplex_volume_intrinsic(d_degen) is None

Computing Dihedral Angles
-------------------------

The dihedral angle is the angle between two faces sharing an edge, measured
from the interior of the dihedral.

.. code-block:: python

   import numpy as np
   from psge.core.geometry_ext import GeometryExtrinsic

   geo = GeometryExtrinsic(dimension=3)

   # Two co-planar faces → dihedral = π (flat geometry)
   p1 = np.array([0., 0., 0.])
   p2 = np.array([1., 0., 0.])
   p3 = np.array([0., 1., 0.])   # above shared edge
   p4 = np.array([0., -1., 0.])  # below shared edge

   angle = geo.dihedral_angle(p1, p2, p3, p4)
   print(f"Flat dihedral = {angle:.4f} rad  (π = {np.pi:.4f})")

   # Two orthogonal faces → dihedral = π/2
   p4_orth = np.array([0., 0., 1.])
   angle_orth = geo.dihedral_angle(p1, p2, p3, p4_orth)
   print(f"Right-angle dihedral = {angle_orth:.4f} rad  (π/2 = {np.pi/2:.4f})")

Computing Regge Deficits and the Regge Action
---------------------------------------------

The *Regge deficit* at an edge measures the discrete curvature:

.. math::

   \varepsilon_h = 2\pi - \sum_{\sigma \supset h} \theta_h(\sigma)

where the sum runs over all simplices :math:`\sigma` containing the hinge
(edge) :math:`h`, and :math:`\theta_h(\sigma)` is the dihedral angle of
:math:`\sigma` at :math:`h`.

.. code-block:: python

   import numpy as np
   from psge.curvature.deficit import deficit_extrinsic
   from psge.curvature.action import regge_action

   # Flat geometry: six equilateral triangles around an interior edge
   angles_flat = np.full(6, np.pi / 3)
   eps = deficit_extrinsic(angles_flat)
   print(f"Flat deficit = {eps:.6f}")   # 0.000000

   # Curved cone: four triangles meeting at the apex
   cone_angle = np.pi * 3 / 8        # each dihedral
   angles_cone = np.full(4, cone_angle)
   eps_cone = deficit_extrinsic(angles_cone)
   print(f"Cone deficit = {eps_cone:.6f}")   # π/2 ≈ 1.5708

   # Regge action  S = Σ deficit_i × dual_volume_i
   deficits = np.array([0.0, eps_cone, 0.0])
   volumes  = np.array([0.5, 0.25, 0.5])
   S = regge_action(deficits, volumes)
   print(f"Regge action S = {S:.6f}")

Using Validation Oracles
------------------------

PSGE-II ships with analytical reference values for testing:

.. code-block:: python

   from psge.validation.oracles import OracleV11

   # Flat geometry must have zero deficit
   assert OracleV11.flat_square_deficit() == 0.0

   # Cone oracle (default π/2 deficit)
   cone_def = OracleV11.cone_deficit()
   print(f"Cone deficit oracle = {cone_def:.4f}")  # π/2 ≈ 1.5708

   # Scaling invariance meta-data
   info = OracleV11.scaling_invariance(scale_factor=3.0)
   print(info)  # {'deficit_invariant': True, 'volume_scales_by': 9.0}
