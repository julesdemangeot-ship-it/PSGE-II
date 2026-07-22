Mathematical Theory
===================

.. contents::
   :local:
   :depth: 2

Regge Calculus
--------------

Regge Calculus (T. Regge, 1961) is a discrete formulation of General
Relativity in which continuous spacetime is replaced by a *piecewise-flat*
simplicial complex.  All curvature is concentrated on co-dimension-2 *hinges*
(edges in 3D, triangles in 4D) as *angular deficits*.

Simplicial Complexes
~~~~~~~~~~~~~~~~~~~~

A *d-simplex* :math:`\sigma^d` is the convex hull of :math:`d+1`
affinely independent points (vertices).  A *simplicial complex* is a
finite collection of simplices closed under taking faces.

In PSGE-II:

* ``0-simplex`` — vertex
* ``1-simplex`` — edge (length :math:`\ell_{ij}`)
* ``2-simplex`` — triangle (face)
* ``3-simplex`` — tetrahedron (cell)

The geometry is encoded entirely in the set of edge lengths
:math:`\{\ell_{ij}\}`, with no reference to an ambient embedding space.

Gram Matrices
~~~~~~~~~~~~~

For a *d*-simplex with vertices :math:`p_0, \ldots, p_d` embedded in
:math:`\mathbb{R}^n`, the **Gram matrix** is

.. math::

   G_{ij} = \langle p_i - p_0,\; p_j - p_0 \rangle, \quad i,j = 0,\ldots,d.

Note that :math:`G_{0j} = G_{i0} = 0` for all :math:`i,j` by construction.
The Gram matrix is symmetric positive-semi-definite for non-degenerate
simplices.

Cayley-Menger Determinants
~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a pairwise distance matrix :math:`D` of :math:`n` points, the
**Cayley-Menger matrix** of size :math:`(n+1)\times(n+1)` is

.. math::

   \mathrm{CM} =
   \begin{pmatrix}
       0 & 1 & 1 & \cdots & 1 \\
       1 & 0 & d_{01}^2 & \cdots & d_{0,n-1}^2 \\
       1 & d_{01}^2 & 0 & \cdots & d_{1,n-1}^2 \\
       \vdots & \vdots & \vdots & \ddots & \vdots \\
       1 & d_{0,n-1}^2 & d_{1,n-1}^2 & \cdots & 0
   \end{pmatrix}.

The volume of the :math:`(n-1)`-simplex spanned by the :math:`n` points
satisfies

.. math::

   \bigl((n-1)!\,V\bigr)^2 = \frac{(-1)^n}{2^{n-1}} \det(\mathrm{CM}).

Dihedral Angles
~~~~~~~~~~~~~~~

The **dihedral angle** :math:`\theta_h(\sigma)` of simplex :math:`\sigma`
at hinge :math:`h` is the angle between the two co-dimension-1 faces of
:math:`\sigma` that share :math:`h`, measured from the interior.

For the **extrinsic** (v1.1) formulation, given outward unit normals
:math:`\hat{n}_1` and :math:`\hat{n}_2` of the two faces:

.. math::

   \theta = \arccos(-\hat{n}_1 \cdot \hat{n}_2).

Key special cases:

* Co-planar faces (flat geometry): :math:`\theta = \pi`.
* Orthogonal faces: :math:`\theta = \pi/2`.

Regge Deficit
~~~~~~~~~~~~~

The **Regge deficit** (angular deficit) at hinge :math:`h` is

.. math::

   \varepsilon_h = 2\pi - \sum_{\sigma \supset h} \theta_h(\sigma),

where the sum is over all simplices meeting at :math:`h`.

* :math:`\varepsilon_h = 0` — locally flat geometry.
* :math:`\varepsilon_h > 0` — positive (elliptic) curvature.
* :math:`\varepsilon_h < 0` — negative (hyperbolic) curvature.

Regge Action
~~~~~~~~~~~~

The discrete Einstein-Hilbert action is

.. math::

   S = \sum_{h} \varepsilon_h \cdot A_h,

where :math:`A_h` is the *dual volume* (area/length) associated with hinge
:math:`h`.  In the continuum limit this converges to the Einstein-Hilbert
action :math:`S = \int R \,\sqrt{g}\,\mathrm{d}^n x`.

v1.1 vs v1.2 Formulations
--------------------------

v1.1 — Extrinsic
~~~~~~~~~~~~~~~~

The v1.1 engine computes all geometric quantities from explicit vertex
coordinates embedded in Euclidean :math:`\mathbb{R}^d`.  This is
straightforward but requires an external ambient space.

v1.2 — Intrinsic (in development)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The v1.2 engine works *purely* from edge lengths, following Regge's original
construction.  Dihedral angles are derived from Cayley-Menger determinants
with no reference to an embedding.  This is the mathematically correct
setting for discrete General Relativity.

References
----------

* T. Regge, *General Relativity without coordinates*, Nuovo Cimento **19**,
  558–571 (1961).
* R. Williams & P. Tuckey, *Regge calculus: a brief review and bibliography*,
  Class. Quantum Grav. **9**, 1409 (1992).
* H. W. Hamber, *Quantum Gravitation*, Springer, 2009.
