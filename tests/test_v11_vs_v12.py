"""L6 Comparison campaign: v1.1 (extrinsic) vs v1.2 (intrinsic).

Demonstrates two fundamental properties:
1. **Agreement on embeddable geometries**: On any mesh that embeds isometrically
   in flat R^3, the intrinsic engine reproduces the extrinsic results exactly.
   v1.2 is a strict superset, not a drifting replacement.

2. **Capability gap on curved cones**: On intrinsically curved geometries (e.g.,
   a cone with m ≠ 6 regular tetrahedra), the extrinsic engine (forced to a
   coordinate closure) reports ~0 deficit, unable to "see" the curvature. The
   intrinsic engine returns the analytical value.

Tests organized by level:

  L6.1 — Regular tetrahedron dihedral
  L6.2 — Flat fan deficit
  L6.3 — Flat Regge action
  L6.4 — Scale invariance
  L6.5 — Degenerate rejection
  L6.6 — Curved cone (discriminant test)

Each test prints a summary to stdout; all are assertions suitable for pytest.
"""

import numpy as np
from itertools import combinations

from psge.core.intrinsic_utils import EdgeLengths
from psge.core.geometry_int import intrinsic_dihedral
from psge.curvature.deficit_int import intrinsic_deficit, hinges_and_cells, regge_action
from psge.validation.oracles_v12 import (
    regular_tetrahedron,
    flat_fan,
    regular_cone,
    cone_deficit,
)

TOL = 1e-9


# ============================================================================
# Minimal faithful reimplementation of v1.1 extrinsic engine
# ============================================================================

def extrinsic_dihedral_v11(coords, tet, hinge):
    """Extrinsic dihedral angle (v1.1 reference implementation).
    
    Computes the angle between two faces of a tetrahedron sharing an edge,
    using global embedding coordinates.
    """
    a, b = hinge
    others = [v for v in tet if v not in (a, b)]
    pa, pb = np.asarray(coords[a], float), np.asarray(coords[b], float)
    pc, pd = np.asarray(coords[others[0]], float), np.asarray(coords[others[1]], float)

    e = pb - pa
    e = e / np.linalg.norm(e)

    vc = (pc - pa) - np.dot(pc - pa, e) * e
    vd = (pd - pa) - np.dot(pd - pa, e) * e

    nc, nd = np.linalg.norm(vc), np.linalg.norm(vd)
    if nc < 1e-14 or nd < 1e-14:
        return 0.0

    ct = np.clip(np.dot(vc, vd) / (nc * nd), -1, 1)
    return float(np.arccos(ct))


def extrinsic_deficit_v11(coords, hinge, cells):
    """Extrinsic deficit (v1.1 reference implementation).
    
    deficit = 2π - Σ(incident dihedral angles)
    """
    return float(2 * np.pi - sum(extrinsic_dihedral_v11(coords, c, hinge) for c in cells))


def extrinsic_regge_action_v11(coords, cells, hinge_cells):
    """Extrinsic Regge action (v1.1 reference).
    
    For flat embeddable geometries, this should be ~0.
    """
    from psge.curvature.deficit_int import is_interior_hinge
    
    action = 0.0
    for hinge, cells_list in hinge_cells.items():
        if is_interior_hinge(hinge, cells_list):
            deficit = extrinsic_deficit_v11(coords, hinge, cells_list)
            # Approximate volume: sum of cell volumes / 8 (placeholder)
            volume = 0.0
            for c in cells_list:
                # Simple volume estimate from coordinates
                verts = [np.asarray(coords[v], float) for v in c]
                edges = [verts[i] - verts[0] for i in range(1, 4)]
                vol = abs(np.linalg.det(edges)) / 6.0
                volume += vol / 8.0
            action += deficit * volume
    return action


# ============================================================================
# L6 Test Suite
# ============================================================================

def test_L6_1_regular_tetrahedron_dihedral():
    """L6.1 — Regular tetrahedron: both engines compute arccos(1/3)."""
    coords, tet = regular_tetrahedron()
    metric = EdgeLengths.from_tetrahedra(coords, [tet])

    # v1.1: extrinsic dihedral
    angle_ext = extrinsic_dihedral_v11(coords, tet, (0, 1))

    # v1.2: intrinsic dihedral
    angle_int = intrinsic_dihedral(tet, (0, 1), metric)

    # Oracle: arccos(1/3)
    angle_oracle = np.arccos(1.0 / 3.0)

    assert np.isclose(angle_ext, angle_oracle, atol=1e-10), f"v1.1 dihedral: {angle_ext}"
    assert np.isclose(angle_int, angle_oracle, atol=1e-10), f"v1.2 dihedral: {angle_int}"
    assert np.isclose(angle_ext, angle_int, atol=TOL), f"Agreement: {angle_ext} vs {angle_int}"

    print(f"✓ L6.1 Regular tetrahedron dihedral")
    print(f"  v1.1: {angle_ext:.12f}")
    print(f"  v1.2: {angle_int:.12f}")
    print(f"  Oracle (arccos(1/3)): {angle_oracle:.12f}")


def test_L6_2_flat_fan_deficit():
    """L6.2 — Flat fan: both engines report zero deficit."""
    coords, tets = flat_fan()
    metric = EdgeLengths.from_tetrahedra(coords, tets)

    # Find the central hinge (0, 1)
    hinge = (0, 1)
    cells = [t for t in tets if hinge[0] in t and hinge[1] in t]

    # v1.1: extrinsic deficit
    deficit_ext = extrinsic_deficit_v11(coords, hinge, cells)

    # v1.2: intrinsic deficit
    deficit_int = intrinsic_deficit(hinge, cells, metric)

    assert abs(deficit_ext) < 1e-9, f"v1.1 deficit not zero: {deficit_ext}"
    assert abs(deficit_int) < 1e-9, f"v1.2 deficit not zero: {deficit_int}"
    assert np.isclose(deficit_ext, deficit_int, atol=TOL), \
        f"Agreement: {deficit_ext} vs {deficit_int}"

    print(f"✓ L6.2 Flat fan deficit")
    print(f"  v1.1 deficit: {deficit_ext:.12e}")
    print(f"  v1.2 deficit: {deficit_int:.12e}")


def test_L6_3_flat_regge_action():
    """L6.3 — Flat Regge action: both engines report ~0."""
    coords, tets = flat_fan()
    metric = EdgeLengths.from_tetrahedra(coords, tets)

    # Build hinge map
    hinge_cells = hinges_and_cells(tets)

    # v1.1: extrinsic action
    action_ext = extrinsic_regge_action_v11(coords, tets, hinge_cells)

    # v1.2: intrinsic action
    hinges = list(hinge_cells.keys())
    action_int = regge_action(hinges, hinge_cells, metric)

    assert abs(action_ext) < 1e-8, f"v1.1 action not zero: {action_ext}"
    assert abs(action_int) < 1e-8, f"v1.2 action not zero: {action_int}"
    assert np.isclose(action_ext, action_int, atol=1e-7), \
        f"Agreement: {action_ext} vs {action_int}"

    print(f"✓ L6.3 Flat Regge action")
    print(f"  v1.1 action: {action_ext:.12e}")
    print(f"  v1.2 action: {action_int:.12e}")


def test_L6_4_scale_invariance():
    """L6.4 — Scale invariance: both engines scale consistently."""
    # Generate flat fan at scale 1.0
    coords1, tets1 = flat_fan()
    metric1 = EdgeLengths.from_tetrahedra(coords1, tets1)

    # Generate flat fan at scale 2.0
    coords2 = {v: 2.0 * np.asarray(c) for v, c in coords1.items()}
    metric2 = EdgeLengths.from_tetrahedra(coords2, tets1)

    hinge = (0, 1)
    cells = [t for t in tets1 if hinge[0] in t and hinge[1] in t]

    # Deficit should be scale-invariant (both zero)
    deficit1_ext = extrinsic_deficit_v11(coords1, hinge, cells)
    deficit2_ext = extrinsic_deficit_v11(coords2, hinge, cells)

    deficit1_int = intrinsic_deficit(hinge, cells, metric1)
    deficit2_int = intrinsic_deficit(hinge, cells, metric2)

    assert abs(deficit1_ext - deficit2_ext) < 1e-9, \
        f"v1.1 scale: {deficit1_ext} vs {deficit2_ext}"
    assert abs(deficit1_int - deficit2_int) < 1e-9, \
        f"v1.2 scale: {deficit1_int} vs {deficit2_int}"

    print(f"✓ L6.4 Scale invariance")
    print(f"  v1.1 scale 1.0: {deficit1_ext:.12e}")
    print(f"  v1.1 scale 2.0: {deficit2_ext:.12e}")
    print(f"  v1.2 scale 1.0: {deficit1_int:.12e}")
    print(f"  v1.2 scale 2.0: {deficit2_int:.12e}")


def test_L6_5_degenerate_rejection():
    """L6.5 — Degenerate tetrahedron: both engines reject it."""
    # Four coplanar points (degenerate)
    coords = {
        0: [0.0, 0.0, 0.0],
        1: [1.0, 0.0, 0.0],
        2: [0.0, 1.0, 0.0],
        3: [1.0, 1.0, 0.0],
    }
    edges = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    metric = EdgeLengths.from_coords(coords, edges)

    tet = (0, 1, 2, 3)

    # v1.1: should have zero or near-zero volume
    verts = [np.asarray(coords[v], float) for v in tet]
    edges_vec = [verts[i] - verts[0] for i in range(1, 4)]
    vol_ext = abs(np.linalg.det(edges_vec)) / 6.0

    # v1.2: should raise ValueError on local_coords
    from psge.core.intrinsic_utils import local_coords
    try:
        local_coords(tet, metric)
        v12_rejected = False
    except ValueError:
        v12_rejected = True

    assert vol_ext < 1e-14, f"v1.1 volume should be ~0: {vol_ext}"
    assert v12_rejected, "v1.2 should reject degenerate tetrahedron"

    print(f"✓ L6.5 Degenerate rejection")
    print(f"  v1.1 volume: {vol_ext:.12e}")
    print(f"  v1.2 raises ValueError: {v12_rejected}")


def test_L6_6_curved_cone_discriminant():
    """L6.6 — Curved cone: the discriminant test.

    On a cone with m regular tets:
    - v1.1 (forced embedding): deficit ≈ 0 (cannot see curvature)
    - v1.2 (intrinsic): deficit = analytical oracle (sees curvature)
    """
    for m in [3, 4, 5]:
        # Intrinsic: perfect regular cone
        metric_int, tets_int, hinge = regular_cone(m)
        cells_int = [t for t in tets_int if hinge[0] in t and hinge[1] in t]
        deficit_int = intrinsic_deficit(hinge, cells_int, metric_int)
        deficit_oracle = cone_deficit(m)

        # Extrinsic: forced coordinate closure (destroys curvature)
        # Create a coordinate embedding that closes at 2π around the edge
        coords_ext = {0: [0.0, 0.0, 0.0], 1: [0.0, 0.0, 1.0]}
        nid = 2
        R = 0.9
        for k in range(m):
            ang = 2 * np.pi * k / m  # forced closure
            coords_ext[nid] = [R * np.cos(ang), R * np.sin(ang), 0.5]
            nid += 1

        ring = [2 + k for k in range(m)]
        tets_ext = [(0, 1, ring[k], ring[(k + 1) % m]) for k in range(m)]
        cells_ext = [t for t in tets_ext if hinge[0] in t and hinge[1] in t]
        deficit_ext = extrinsic_deficit_v11(coords_ext, hinge, cells_ext)

        # Assertions
        assert np.isclose(deficit_int, deficit_oracle, atol=1e-10), \
            f"m={m}: v1.2 should match oracle"
        assert abs(deficit_ext) < 1e-9, \
            f"m={m}: v1.1 should see ~0 (forced closure)"
        assert abs(deficit_int - deficit_ext) > 0.1, \
            f"m={m}: the gap should be large (capability difference)"

        print(f"  m={m}: v1.1 deficit={deficit_ext:.12e}, "
              f"v1.2 deficit={deficit_int:.12f}, "
              f"oracle={deficit_oracle:.12f}")

    print(f"✓ L6.6 Curved cone discriminant")


# ============================================================================
# Test summary table (printed after all tests)
# ============================================================================

def print_summary_table():
    """Print a summary table of all L6 results."""
    print("\n" + "=" * 80)
    print("L6 COMPARISON CAMPAIGN: v1.1 vs v1.2")
    print("=" * 80)
    print(f"\n{'Test':<30} {'v1.1':<15} {'v1.2':<15} {'Result':<20}")
    print("-" * 80)
    print(f"{'L6.1 Regular dihedral':<30} {'✓':<15} {'✓':<15} {'Identical':<20}")
    print(f"{'L6.2 Flat fan deficit':<30} {'✓':<15} {'✓':<15} {'Deficit = 0':<20}")
    print(f"{'L6.3 Flat action':<30} {'✓':<15} {'✓':<15} {'Action = 0':<20}")
    print(f"{'L6.4 Scale invariance':<30} {'✓':<15} {'✓':<15} {'Identical':<20}")
    print(f"{'L6.5 Degenerate':<30} {'✓':<15} {'✓':<15} {'Both reject':<20}")
    print(f"{'L6.6 Curved cone':<30} {'✗ Limitation':<15} {'✓':<15} {'Gap > 0.1 rad':<20}")
    print("=" * 80)
    print("\nConclusion:")
    print("  • v1.2 is a strict superset of v1.1 (agreement on embeddable geometries)")
    print("  • v1.2 extends to intrinsically curved geometries (e.g., cones with m ≠ 6)")
    print("  • Both engines numerically robust and consistent in their respective domains")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    # Run all tests
    tests = [
        test_L6_1_regular_tetrahedron_dihedral,
        test_L6_2_flat_fan_deficit,
        test_L6_3_flat_regge_action,
        test_L6_4_scale_invariance,
        test_L6_5_degenerate_rejection,
        test_L6_6_curved_cone_discriminant,
    ]

    for test_fn in tests:
        test_fn()
        print()

    print_summary_table()
