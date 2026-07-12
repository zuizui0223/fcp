"""General multi-patch extension of the FCP phase-boundary model.

For n patches, let p_i be the frequency of morph A. Local directional
selection is d + h_i, within-patch balancing strength is b, and migration
occurs on an undirected weighted graph with Laplacian L at overall scale m.

Near the all-B boundary, rare A grows according to

    x' = [(b + d) I + diag(h) - m L] x.

Near the all-A boundary, rare B grows according to

    y' = [(b - d) I + diag(-h) - m L] y.

Define

    H_A = lambda_max(diag(h)  - m L)
    H_B = lambda_max(diag(-h) - m L)

and the spectral compression

    H_bal = (H_A + H_B) / 2
    D_asym = (H_A - H_B) / 2.

Then

    lambda_A|B = b + H_bal + d + D_asym
    lambda_B|A = b + H_bal - (d + D_asym)

so protected polymorphism exists iff

    b + H_bal > abs(d + D_asym).

For the symmetric two-patch model h=(+q,-q), this reduces exactly to

    H_bal = sqrt(q^2 + m^2) - m,  D_asym = 0.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Iterable, List, Sequence, Tuple


Matrix = List[List[float]]


@dataclass(frozen=True)
class MultiPatchParameters:
    b: float
    d: float
    h: Tuple[float, ...]
    adjacency: Tuple[Tuple[float, ...], ...]
    m: float

    def __post_init__(self) -> None:
        n = len(self.h)
        if n < 1:
            raise ValueError("at least one patch is required")
        if self.m < 0:
            raise ValueError("m must be non-negative")
        if len(self.adjacency) != n or any(len(row) != n for row in self.adjacency):
            raise ValueError("adjacency must be square and match h")
        for i in range(n):
            if abs(self.adjacency[i][i]) > 1e-12:
                raise ValueError("adjacency diagonal must be zero")
            for j in range(n):
                if self.adjacency[i][j] < 0:
                    raise ValueError("adjacency weights must be non-negative")
                if abs(self.adjacency[i][j] - self.adjacency[j][i]) > 1e-12:
                    raise ValueError("adjacency must be symmetric")


def graph_laplacian(adjacency: Sequence[Sequence[float]]) -> Matrix:
    """Return the weighted graph Laplacian D - A."""
    n = len(adjacency)
    if any(len(row) != n for row in adjacency):
        raise ValueError("adjacency must be square")
    L = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        degree = 0.0
        for j in range(n):
            w = float(adjacency[i][j])
            if w < 0:
                raise ValueError("adjacency weights must be non-negative")
            degree += w
            if i != j:
                L[i][j] = -w
        L[i][i] = degree
    return L


def selection_migration_matrix(
    h: Sequence[float], adjacency: Sequence[Sequence[float]], m: float
) -> Matrix:
    """Return diag(h) - m L for an undirected migration network."""
    if m < 0:
        raise ValueError("m must be non-negative")
    n = len(h)
    if len(adjacency) != n or any(len(row) != n for row in adjacency):
        raise ValueError("adjacency must be square and match h")
    L = graph_laplacian(adjacency)
    return [
        [((float(h[i]) if i == j else 0.0) - m * L[i][j]) for j in range(n)]
        for i in range(n)
    ]


def _largest_eigenvalue_symmetric(matrix: Sequence[Sequence[float]], tol: float = 1e-12) -> float:
    """Largest eigenvalue of a real symmetric matrix via Jacobi rotations.

    This dependency-free routine is intended for modest research-prototype
    matrices. It validates symmetry and diagonalizes a copy of the matrix.
    """
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("matrix must be non-empty and square")
    a = [[float(x) for x in row] for row in matrix]
    for i in range(n):
        for j in range(n):
            if abs(a[i][j] - a[j][i]) > 1e-10:
                raise ValueError("matrix must be symmetric")
    if n == 1:
        return a[0][0]

    max_iter = 100 * n * n
    for _ in range(max_iter):
        p, q = 0, 1
        off = abs(a[p][q])
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > off:
                    p, q = i, j
                    off = abs(a[i][j])
        if off < tol:
            break

        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        tau = (aqq - app) / (2.0 * apq)
        t = (1.0 if tau >= 0 else -1.0) / (abs(tau) + sqrt(1.0 + tau * tau))
        c = 1.0 / sqrt(1.0 + t * t)
        s = t * c

        for k in range(n):
            if k == p or k == q:
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = c * akp - s * akq
            a[k][q] = a[q][k] = s * akp + c * akq

        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        a[p][q] = a[q][p] = 0.0
    else:
        raise RuntimeError("Jacobi eigensolver did not converge")

    return max(a[i][i] for i in range(n))


def spectral_components(params: MultiPatchParameters) -> Tuple[float, float, float, float]:
    """Return (H_A, H_B, H_bal, D_asym)."""
    matrix_a = selection_migration_matrix(params.h, params.adjacency, params.m)
    matrix_b = selection_migration_matrix(tuple(-x for x in params.h), params.adjacency, params.m)
    h_a = _largest_eigenvalue_symmetric(matrix_a)
    h_b = _largest_eigenvalue_symmetric(matrix_b)
    h_bal = 0.5 * (h_a + h_b)
    d_asym = 0.5 * (h_a - h_b)
    return h_a, h_b, h_bal, d_asym


def invasion_rates(params: MultiPatchParameters) -> Tuple[float, float]:
    """Return dominant rare-invasion exponents (A into B, B into A)."""
    h_a, h_b, _, _ = spectral_components(params)
    return params.b + params.d + h_a, params.b - params.d + h_b


def classify_phase(params: MultiPatchParameters, tol: float = 1e-12) -> str:
    """Classify boundary invasibility into four qualitative phases."""
    a_into_b, b_into_a = invasion_rates(params)
    a_pos = a_into_b > tol
    b_pos = b_into_a > tol
    if a_pos and b_pos:
        return "protected_polymorphism"
    if a_pos and not b_pos:
        return "A_fixation"
    if not a_pos and b_pos:
        return "B_fixation"
    return "bistability"


def complete_graph(n: int, weight: float = 1.0) -> Tuple[Tuple[float, ...], ...]:
    """Convenience constructor for an undirected complete graph."""
    if n < 1 or weight < 0:
        raise ValueError("n must be positive and weight non-negative")
    return tuple(
        tuple(0.0 if i == j else float(weight) for j in range(n))
        for i in range(n)
    )


def path_graph(n: int, weight: float = 1.0) -> Tuple[Tuple[float, ...], ...]:
    """Convenience constructor for an undirected path graph."""
    if n < 1 or weight < 0:
        raise ValueError("n must be positive and weight non-negative")
    rows = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n - 1):
        rows[i][i + 1] = rows[i + 1][i] = float(weight)
    return tuple(tuple(row) for row in rows)
