"""Numerical realizability checks for certified Gram-equation solutions."""

from __future__ import annotations

from math import cos, pi
from typing import List, Sequence, Tuple


def gram_matrix(level: int, numerators: Sequence[int]) -> List[List[float]]:
    if len(numerators) != 6:
        raise ValueError("six angle numerators are required")
    # Outward face normals form angle pi-alpha across an edge, hence their
    # scalar product is cos(pi-alpha)=-cos(alpha).
    c12, c34, c13, c24, c14, c23 = (
        -cos(pi * m / level) for m in numerators
    )
    return [
        [1.0, c12, c13, c14],
        [c12, 1.0, c23, c24],
        [c13, c23, 1.0, c34],
        [c14, c24, c34, 1.0],
    ]


def _det3(a: Sequence[Sequence[float]]) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def _minor3(g: Sequence[Sequence[float]], skip_row: int, skip_col: int) -> float:
    rows = [i for i in range(4) if i != skip_row]
    cols = [j for j in range(4) if j != skip_col]
    return _det3([[g[i][j] for j in cols] for i in rows])


def cofactors(g: Sequence[Sequence[float]]) -> List[List[float]]:
    return [
        [((-1.0) ** (i + j)) * _minor3(g, i, j) for j in range(4)]
        for i in range(4)
    ]


def gram_determinant_from_cosines(c: Sequence[float]) -> float:
    """Expanded determinant in paper edge order.

    This form also exposes the quadratic dependence on the last cosine.
    """

    a, b, c13, d, e, f = c  # c12,c34,c13,c24,c14,c23
    return (
        1.0
        - (a * a + b * b + c13 * c13 + d * d + e * e + f * f)
        + 2.0 * (a * c13 * f + a * e * d + c13 * e * b + f * d * b)
        + a * a * b * b
        + c13 * c13 * d * d
        + e * e * f * f
        - 2.0 * (a * c13 * d * b + a * e * f * b + c13 * e * f * d)
    )


def is_nondegenerate_tetrahedron(
    level: int,
    numerators: Sequence[int],
    tolerance: float = 1e-10,
) -> bool:
    """Test the paper's type-(iv) condition for a certified solution.

    All four diagonal cofactors must be positive (every triple of normals is
    independent), and the unique null relation must have coefficients of one
    sign.  The caller should separately certify det(G)=0 exactly.
    """

    g = gram_matrix(level, numerators)
    c = cofactors(g)
    if not all(c[i][i] > tolerance for i in range(4)):
        return False
    # Column 0 of adj(G), equivalently row 0 of the cofactor matrix, is a
    # null vector.  Its first entry is positive by the preceding test.
    null_relation = [c[0][j] for j in range(4)]
    return all(value > tolerance for value in null_relation)


def quadratic_coefficients_for_c23(c: Sequence[float]) -> Tuple[float, float, float]:
    """Return A,B,C with det(G)=A*c23^2+B*c23+C.

    Input order is c12,c34,c13,c24,c14.
    """

    a, b, c13, d, e = c
    A = e * e - 1.0
    B = 2.0 * (a * c13 + d * b - a * e * b - c13 * e * d)
    C = (
        1.0
        - (a * a + b * b + c13 * c13 + d * d + e * e)
        + 2.0 * (a * e * d + c13 * e * b)
        + a * a * b * b
        + c13 * c13 * d * d
        - 2.0 * a * c13 * d * b
    )
    return A, B, C
