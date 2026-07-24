# PSGE-II

Intrinsic Euclidean Regge Geometry Engine

## Overview

PSGE-II is a computational geometry engine implementing the Regge Calculus formalism for discrete spacetime geometry. The project is organized in progressive versions, each adding capabilities while maintaining analytical validation.

## Versions

### PSGE-II v1.1 (Stable Reference)

**Status:** ✅ Stable and validated

Extrinsic Euclidean engine validated through comprehensive analytical oracle campaigns:

- **Dihedral angle computation** from global coordinates
- **Regge deficits** on embeddable meshes (flat geometry)
- **Regge action** in flat Euclidean geometry
- **Geometric invariances:** isometry, scale invariance
- **Degenerate case detection** and handling

This version serves as the stable reference for the current engine.

### PSGE-II v1.2 (Recommended Evolution)

**Status:** 🚀 In development on branch `develop/v1.2`

Transition to an intrinsic formulation that removes the dependency on external embeddings:

- **Intrinsic dihedral angles:** computed from edge lengths alone via Gram/Cayley-Menger determinants
- **Intrinsic Regge deficits:** derived purely from metric data
- **Barycentric dual measures:** validation of dual geometry computations
- **Curved cone oracle:** introduction of cone geometries with non-zero deficit
- **Comprehensive validation campaign:** full suite of tests for intrinsically curved geometries

This version extends v1.1 capabilities and aligns with the classical Regge Calculus construction.

### PSGE-II v2-dev (Research Branch)

**Status:** 🔬 Research and exploration

Generalization to Lorentzian metrics and spacetime geometry:

- **Lorentzian dihedral angles:** angles and rapidities based on face signature
- **Lorentzian Gram matrices:** signature-preserving metric formulation
- **Signature-adapted dual measures:** generalized barycentric constructs
- **Lorentzian Regge action:** extension of action principle to spacetime
- **Dedicated validation campaign:** thorough testing before stable promotion

## Repository Structure

```
PSGE-II/
├── src/
│   └── psge/
│       ├── core/              # Geometric and mathematical kernel
│       ├── mesh/              # Mesh management and topology
│       ├── curvature/         # Curvature analysis and actions
│       └── validation/        # Test campaigns and oracles
├── tests/                     # Unit tests (pytest)
├── docs/                      # Documentation and theory
├── reports/                   # Generated validation reports
└── pyproject.toml            # Package metadata
```

## Quick Start

### Installation

```bash
git clone https://github.com/julesdemangeot-ship-it/PSGE-II.git
cd PSGE-II
pip install -e .
```

### Running Tests

```bash
pytest tests/
```

### Generating Validation Reports

```bash
python -m psge.validation.suite
```

## Development

- **Stable releases:** `stable/v1.1` (protected branch)
- **Active development:** `develop/v1.2` for intrinsic formulation
- **Main branch:** integration point for validated features

## Contributing

See `CONTRIBUTING.md` for guidelines on contributing to the project.

## License

See `LICENSE` file.

## References

- Regge, T. (1961). "General Relativity Without Coordinates"
- Williams, R.M. & Tuckey, P.A. (1992). "Regge Calculus: A Brief Review and Bibliography"

## Don

Soutenez mon travail ici : [Ko-fi](https://ko-fi.com/julesdemangeot)
