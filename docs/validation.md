# Validation Strategy

## v1.1 Validation Campaign

The v1.1 extrinsic engine has been validated against analytical oracles covering:

### Test Categories

1. **Dihedral Angle Computation**
   - Flat configurations (π expected)
   - Acute/obtuse configurations
   - Degenerate cases

2. **Regge Deficit**
   - Flat meshes (deficit = 0)
   - Cone geometries (non-zero deficit)
   - Validation with analytical reference values

3. **Regge Action**
   - Flat spacetime (action = 0)
   - Curved configurations
   - Consistency checks

4. **Geometric Invariances**
   - Isometry preservation
   - Scale invariance
   - Reparameterization invariance

5. **Degenerate Case Detection**
   - Coplanar configurations
   - Zero-volume simplices
   - Ill-conditioned metrics

## v1.2 Validation Plan

### New Test Categories

1. **Intrinsic Dihedral Angles**
   - Comparison with v1.1 extrinsic computation
   - Cayley-Menger determinant validation
   - Gram matrix consistency

2. **Intrinsic Deficit Computation**
   - Validation against embedded reference
   - Handling of non-embeddable configurations

3. **Curved Cone Oracles**
   - Cone geometries with specified deficit
   - Non-zero deficit validation
   - Curvature singularities

4. **Barycentric Dual Measures**
   - Dual volume computation
   - Topological consistency
   - Sum rules verification

## Validation Infrastructure

### Oracle System
Analytical reference values for key test cases:
- Flat simplex: deficit = 0
- Regular tetrahedron: known deficit values
- Cone configurations: computed analytically

### Automated Test Suite
- Pytest-based infrastructure
- Continuous integration on each commit
- Regression testing

### Report Generation
- PDF reports documenting all test results
- Comparison tables (numerical vs analytical)
- Visualization of validated geometries
