# PSGE-II Theory

## Regge Calculus

Regge Calculus is a formulation of General Relativity on a discrete simplicial complex.
It represents spacetime as a triangulated manifold where geometry is encoded in edge lengths.

### Key Concepts

#### Simplicial Complexes
- A simplicial complex is a collection of simplices (points, edges, triangles, tetrahedra)
- Connectivity is defined topologically

#### Dihedral Angles
- The dihedral angle is the angle between two faces sharing an edge
- In extrinsic formulation: computed from embedded coordinates
- In intrinsic formulation: derived purely from edge lengths

#### Regge Deficit
- The deficit at an edge: `deficit = 2π - Σ(angles around edge)`
- For flat geometry: deficit = 0
- For curved geometry: deficit ≠ 0 (measures curvature)

#### Regge Action
- `S = Σ(deficit × dual_volume)` over all edges
- Discrete analog of Einstein-Hilbert action

## v1.1: Extrinsic Formulation

### Overview
v1.1 computes geometry from explicit coordinates embedded in Euclidean space.

### Advantages
- Direct geometric interpretation
- Simple numerical implementation
- Easy to visualize

### Limitations
- Requires external embedding coordinates
- Cannot represent intrinsically curved geometries
- Limited to embeddable configurations

## v1.2: Intrinsic Formulation

### Overview
v1.2 computes geometry purely from edge lengths, following classical Regge construction.

### Methods

#### Gram Matrices
For a simplex with vertices and edge lengths, the Gram matrix encodes all distance information.

#### Cayley-Menger Determinants
Used to compute volumes and geometric properties from distances.

### Advantages
- Fully intrinsic: no external embedding needed
- Can represent intrinsically curved geometries
- Aligns with classical Regge construction

## v2-dev: Lorentzian Extension

### Overview
Extension to Lorentzian metrics for spacetime geometry.

### Key Differences
- Metric signature (1, 3) instead of (4, 0)
- Dihedral angles become rapidities for timelike edges
- Action principle extends to spacetime manifolds
