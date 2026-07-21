"""PSGE-II v1.2 -- Intrinsic Euclidean Regge geometry.

The v1.1 engine computes dihedral angles from *global* vertex coordinates.
That forces every embeddable mesh to be flat at interior hinges: the dihedral
angles around any properly tiled interior edge sum to exactly 2*pi, so the
Regge deficit is identically 0. Genuinely curved simplicial geometries have no
isometric embedding in flat R^3 and are therefore unreachable.

v1.2 removes the embedding. The only geometric input is the intrinsic metric,
i.e. an edge-length function l: (i, j) -> length. Each tetrahedron is embedded
in *its own* local frame reconstructed from *its own* six edge lengths (Gram
matrix -> Cholesky). The dihedral angle at a hinge is then read off in that
local frame. Because every cell uses an independent frame, the dihedral angles
around a shared edge need not sum to 2*pi, so interior deficits can be non-zero
-- which is exactly the classical Regge construction.

References:
    Regge, T. (1961). "General Relativity Without Coordinates".
    Cayley-Menger determinant for simplex volumes from edge lengths.
"""

from __future__ import annotations

import numpy as np
from itertools import combinations
from typing import Dict, Iterable, Sequence, Tuple

Vertex = int
Edge = Tuple[Vertex, Vertex]
Tet = Tuple[Vertex, Vertex, Vertex, Vertex]

_DEGENERATE_TOL = 1e-12


class EdgeLengths:
    """Intrinsic metric of a simplicial complex.

    Holds a symmetric edge-length function. This is the *only* geometric data
    the intrinsic engine consumes -- no global embedding is ever required.
    """

    def __init__(self, lengths: Dict[Edge, float]):
        self._l: Dict[frozenset, float] = {}
        for (i, j), value in lengths.items():
            if i == j:
                raise ValueError(f"self-edge ({i}, {j}) is not a valid edge")
            key = frozenset((i, j))
            v = float(value)
            if v <= 0.0:
                raise ValueError(f"edge {tuple(sorted((i, j)))} has non-positive length {v}")
            self._l[key] = v

    def length(self, i: Vertex, j: Vertex) -> float:
        try:
            return self._l[frozenset((i, j))]
        except KeyError:
            raise KeyError(f"no length defined for edge {tuple(sorted((i, j)))}")

    def sq(self, i: Vertex, j: Vertex) -> float:
        L = self.length(i, j)
        return L * L

    # --- constructors -----------------------------------------------------

    @classmethod
    def from_coords(cls, coords: Dict[Vertex, Sequence[float]],
                    edges: Iterable[Edge]) -> "EdgeLengths":
        """Derive an intrinsic metric by measuring an existing embedding.

        Useful for regression: an embeddable mesh fed through the intrinsic
        engine must reproduce the extrinsic (v1.1) result exactly.
        """
        d: Dict[Edge, float] = {}
        for (i, j) in edges:
            pi = np.asarray(coords[i], dtype=float)
            pj = np.asarray(coords[j], dtype=float)
            d[(i, j)] = float(np.linalg.norm(pi - pj))
        return cls(d)

    @classmethod
    def from_tetrahedra(cls, coords: Dict[Vertex, Sequence[float]],
                        tets: Iterable[Tet]) -> "EdgeLengths":
        edges = set()
        for t in tets:
            for e in combinations(t, 2):
                edges.add(tuple(sorted(e)))
        return cls.from_coords(coords, edges)


# --------------------------------------------------------------------------
# Single-cell intrinsic geometry
# --------------------------------------------------------------------------

def gram_matrix(tet: Tet, base: Vertex, metric: EdgeLengths):
    """3x3 Gram matrix of the edge vectors emanating from `base`.

    Built purely from edge lengths via the polarization identity
        <u, v> = (|u|^2 + |v|^2 - |u - v|^2) / 2.
    Returns (G, others) where `others` are the three non-base vertices in the
    row/column order used by G.
    """
    others = [v for v in tet if v != base]
    if len(others) != 3:
        raise ValueError("gram_matrix expects a tetrahedron (4 distinct vertices)")
    G = np.empty((3, 3))
    for a in range(3):
        G[a, a] = metric.sq(base, others[a])
    for a, b in combinations(range(3), 2):
        val = 0.5 * (metric.sq(base, others[a])
                     + metric.sq(base, others[b])
                     - metric.sq(others[a], others[b]))
        G[a, b] = G[b, a] = val
    return G, others


def local_coords(tet: Tet, metric: EdgeLengths, base: Vertex | None = None
                 ) -> Dict[Vertex, np.ndarray]:
    """Embed a single tetrahedron in its own local frame from edge lengths.

    A single valid tetrahedron always embeds in R^3, so this is well defined
    per cell even when the *global* complex does not embed. Raises ValueError
    on a degenerate (non positive-definite Gram) tetrahedron.
    """
    base = tet[0] if base is None else base
    G, others = gram_matrix(tet, base, metric)
    G = 0.5 * (G + G.T)  # symmetrize for numerical safety
    if np.linalg.det(G) <= _DEGENERATE_TOL * max(1.0, np.trace(G) ** 3):
        raise ValueError("degenerate tetrahedron (vanishing Gram determinant)")
    try:
        Lc = np.linalg.cholesky(G)  # G = Lc @ Lc.T, Lc lower-triangular
    except np.linalg.LinAlgError:
        raise ValueError("degenerate tetrahedron (Gram not positive-definite)")
    coords: Dict[Vertex, np.ndarray] = {base: np.zeros(3)}
    for k, v in enumerate(others):
        coords[v] = Lc[k].copy()
    return coords


def cayley_menger_volume(tet: Tet, metric: EdgeLengths) -> float:
    """Tetrahedron volume from edge lengths via the Cayley-Menger determinant.

    Independent of `local_coords`, so it doubles as a cross-check on the Gram
    construction. Returns 0.0 for degenerate cells.
    """
    idx = list(tet)
    n = 4
    B = np.ones((n + 1, n + 1))
    B[0, 0] = 0.0
    for a in range(n):
        for b in range(n):
            B[a + 1, b + 1] = 0.0 if a == b else metric.sq(idx[a], idx[b])
    det = float(np.linalg.det(B))
    # V^2 = det(B) / (288)  for a tetrahedron (n = 3)
    v2 = det / 288.0
    return float(np.sqrt(v2)) if v2 > 0.0 else 0.0


def _dihedral_from_coords(coords: Dict[Vertex, np.ndarray], hinge: Edge) -> float:
    """Angle between the two faces sharing `hinge`, in [0, pi]."""
    a, b = hinge
    others = [v for v in coords if v not in (a, b)]
    pa, pb = coords[a], coords[b]
    pc, pd = coords[others[0]], coords[others[1]]
    e = pb - pa
    ne = np.linalg.norm(e)
    if ne < 1e-15:
        return 0.0
    e = e / ne
    vc = (pc - pa) - np.dot(pc - pa, e) * e
    vd = (pd - pa) - np.dot(pd - pa, e) * e
    nc, nd = np.linalg.norm(vc), np.linalg.norm(vd)
    if nc < 1e-14 or nd < 1e-14:
        return 0.0
    ct = np.clip(np.dot(vc, vd) / (nc * nd), -1.0, 1.0)
    return float(np.arccos(ct))


def intrinsic_dihedral(tet: Tet, hinge: Edge, metric: EdgeLengths) -> float:
    """Dihedral angle at `hinge` inside `tet`, from edge lengths alone."""
    a, b = hinge
    if a not in tet or b not in tet:
        raise ValueError(f"hinge {hinge} is not an edge of tet {tet}")
    coords = local_coords(tet, metric)
    return _dihedral_from_coords(coords, hinge)


# --------------------------------------------------------------------------
# Curvature
# --------------------------------------------------------------------------

def intrinsic_deficit(hinge: Edge, cells: Sequence[Tet], metric: EdgeLengths) -> float:
    """Regge deficit at `hinge`: 2*pi minus the sum of incident dihedrals.

    Unlike the extrinsic engine this can be non-zero at an interior hinge.
    """
    total = sum(intrinsic_dihedral(c, hinge, metric) for c in cells)
    return float(2.0 * np.pi - total)


# --------------------------------------------------------------------------
# Analytical oracles for the curved regime
# --------------------------------------------------------------------------

def cone_deficit(m: int) -> float:
    """Deficit of m regular tetrahedra glued around a common edge."""
    return float(2.0 * np.pi - m * np.arccos(1.0 / 3.0))


def regular_cone(m: int):
    """Closed loop of m regular unit tetrahedra sharing a central edge (0, 1).

    Purely intrinsic: all edges have length 1, so every cell is a regular
    tetrahedron. For m != 6 this configuration has *no* isometric embedding in
    flat R^3 -- it is exactly the kind of curved geometry the extrinsic engine
    cannot represent. The central edge (0, 1) is topologically interior.

    Returns (metric, tetrahedra, central_hinge).
    """
    if m < 3:
        raise ValueError("a cone needs at least 3 tetrahedra")
    A, B = 0, 1
    ring = [2 + k for k in range(m)]
    tets = [(A, B, ring[k], ring[(k + 1) % m]) for k in range(m)]
    edges = set()
    for t in tets:
        for e in combinations(t, 2):
            edges.add(tuple(sorted(e)))
    metric = EdgeLengths({e: 1.0 for e in edges})
    return metric, tets, (A, B)


# --------------------------------------------------------------------------
# Topology helper (interior vs boundary) -- unchanged logic from v1.1
# --------------------------------------------------------------------------

def is_interior_hinge(hinge: Edge, cells: Sequence[Tet]) -> bool:
    """A hinge is interior iff every incident face is shared by two cells."""
    face_count: Dict[Tuple[int, int, int], int] = {}
    for c in cells:
        for f in combinations(sorted(c), 3):
            face_count[f] = face_count.get(f, 0) + 1
    a, b = hinge
    for c in cells:
        for f in combinations(sorted(c), 3):
            if a in f and b in f and face_count[tuple(sorted(f))] == 1:
                return False
    return True
