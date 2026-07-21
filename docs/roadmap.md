# PSGE-II Development Roadmap

## Current Status

### v1.1 - Stable (Released)
✅ **Status: STABLE AND VALIDATED**

- [x] Extrinsic Euclidean geometry engine
- [x] Dihedral angle computation
- [x] Regge deficit calculation
- [x] Regge action computation
- [x] Geometric invariance checks
- [x] Degenerate case detection
- [x] Comprehensive oracle validation
- [x] Test suite (pytest)
- [x] Documentation

**Branch:** `stable/v1.1` (protected)

### v1.2 - In Development (Current Focus)
🚀 **Status: ACTIVE DEVELOPMENT**

#### Phase 1: Architecture Setup
- [x] Repository structure
- [x] Module organization
- [x] Branch creation (`develop/v1.2`)
- [ ] Development environment configuration

#### Phase 2: Core Intrinsic Implementations
- [ ] Gram matrix computations from edge lengths
- [ ] Cayley-Menger determinant calculations
- [ ] Intrinsic dihedral angle computation
- [ ] Intrinsic Regge deficit computation
- [ ] Numerical stability analysis

#### Phase 3: Extended Features
- [ ] Barycentric dual measure validation
- [ ] Curved cone oracle implementation
- [ ] Non-zero deficit geometry support
- [ ] Degenerate configuration handling

#### Phase 4: Validation and Testing
- [ ] Comparison against v1.1 (on embeddable configs)
- [ ] Curved geometry validation
- [ ] Oracle-based test suite
- [ ] Regression testing
- [ ] Performance benchmarking

#### Phase 5: Release
- [ ] Code freeze and stabilization
- [ ] Final testing campaign
- [ ] Documentation finalization
- [ ] Promotion to stable version

**Branch:** `develop/v1.2` (active development)

**Timeline:** Q3-Q4 2026

### v2-dev - Research (Future)
🔬 **Status: PLANNING**

#### Goals
- Lorentzian metric extension
- Spacetime geometry support
- Dihedral angles and rapidities
- Lorentzian Regge action
- Validation for relativistic configurations

**Branch:** `research/v2-lorentzian` (planned)

**Estimated Start:** Q1 2027

## Branching Strategy

```
main (integration branch)
 ├── stable/v1.1 (protected, reference)
 ├── develop/v1.2 (active, feature development)
 ├── feature/* (as needed for specific tasks)
 └── research/v2-lorentzian (future)
```

## Key Milestones

1. **v1.2 Alpha** (End of August 2026)
   - Core intrinsic computations working
   - Initial test suite

2. **v1.2 Beta** (End of September 2026)
   - All features implemented
   - Validation campaign complete

3. **v1.2 Release** (End of October 2026)
   - Merge to main
   - Tag v1.2.0
   - Update stable reference

4. **v2-dev Initiation** (Q1 2027)
   - Research branch created
   - Lorentzian extension begins
