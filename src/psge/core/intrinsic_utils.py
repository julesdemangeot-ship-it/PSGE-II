"""PSGE-II v1.2 intrinsic geometry utilities.

Low-level helpers for intrinsic metric computation: edge-length storage,
Gram matrices, Cholesky decomposition, local coordinate reconstruction,
and dihedral angle extraction.

These functions are the mathematical backbone; the public API wraps them
in `geometry_int.py`.
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
    """Intrinsic metric: symmetric edge-length function.

    This is the *only* geometric data the intrinsic engine consumes.
    No global embedding is ever required.
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
        """Get edge length (symmetric)."""
        try:
            return self._l[frozenset((i, j))]
        except KeyError:
            raise KeyError(f"no length defined for edge {tuple(sorted((i, j)))}")

    def sq(self, i: Vertex, j: Vertex) -> float:
        """Get squared edge length."""
        L = self.length(i, j)
        return L * L

    @classmethod
    def from_coords(cls, coords: Dict[Vertex, Sequence[float]],
                    edges: Iterable[Edge]) -> "EdgeLengths":
        """Derive intrinsic metric by measuring an embedded configuration.

        Useful for regression testing: embeddable meshes fed through the
        intrinsic engine should reproduce extrinsic results exactly.
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
        """Derive intrinsic metric from tetrahedra and their embedding."""
        edges = set()
        for t in tets:
            for e in combinations(t, 2):
                edges.add(tuple(sorted(e)))
        return cls.from_coords(coords, edges)


def gram_matrix(tet: Tet, base: Vertex, metric: EdgeLengths) -> Tuple[np.ndarray, list]:
    """3x3 Gram matrix of edge vectors from `base` vertex.

    Built purely from edge lengths via the polarization identity:
        <u, v> = (|u|^2 + |v|^2 - |u - v|^2) / 2

    Returns (G, others) where `others` are the three non-base vertices
    in the order used as row/column indices of G.
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

    A valid tetrahedron always embeds in R^3 (even if the global complex
    doesn't), so this is well-defined per cell. Raises ValueError on
    degenerate (non positive-definite Gram) tetrahedra.
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


def dihedral_from_coords(coords: Dict[Vertex, np.ndarray], hinge: Edge) -> float:
    """Dihedral angle between two faces sharing `hinge`, in [0, π].

    Given local coordinates (e.g., from Cholesky decomposition of Gram matrix),
    extract the dihedral angle at the hinge by computing normal vectors
    to the two incident faces.
    """
    a, b = hinge
    others = [v for v in coords if v not in (a, b)]
    pa, pb = coords[a], coords[b]
    pc, pd = coords[others[0]], coords[others[1]]

    e = pb - pa
    ne = np.linalg.norm(e)
    if ne < 1e-15:
        return 0.0

    e = e / ne

    # Project pc and pd onto the plane perpendicular to edge
    vc = (pc - pa) - np.dot(pc - pa, e) * e
    vd = (pd - pa) - np.dot(pd - pa, e) * e

    nc, nd = np.linalg.norm(vc), np.linalg.norm(vd)
    if nc < 1e-14 or nd < 1e-14:
        return 0.0

    ct = np.clip(np.dot(vc, vd) / (nc * nd), -1.0, 1.0)
    return float(np.arccos(ct))


def cayley_menger_volume(tet: Tet, metric: EdgeLengths) -> float:
    """Tetrahedron volume from edge lengths via Cayley-Menger determinant.

    Independent of `local_coords`, so it serves as a cross-check on Gram
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
    v2 = det / 288.0  # V^2 = det(B) / 288 for tetrahedron

    return float(np.sqrt(v2)) if v2 > 0.0 else 0.0
