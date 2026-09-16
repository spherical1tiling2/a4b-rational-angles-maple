"""Exploratory low-order enumeration followed by exact certification."""

from __future__ import annotations

from math import acos, copysign, cos, isfinite, pi, sqrt
from time import perf_counter
from typing import Dict, List, Sequence, Set, Tuple

from .algebra import exact_torsion_zero
from .geometry import (
    gram_determinant_from_cosines,
    is_nondegenerate_tetrahedron,
    quadratic_coefficients_for_c23,
)


def _quadratic_roots(a: float, b: float, c: float) -> Tuple[float, ...]:
    discriminant = b * b - 4.0 * a * c
    if discriminant < -1e-13:
        return ()
    discriminant = max(0.0, discriminant)
    root_d = sqrt(discriminant)
    if abs(a) < 1e-16:
        return () if abs(b) < 1e-16 else (-c / b,)
    # Stable quadratic formula; the second root is recovered from the product.
    q = -0.5 * (b + copysign(root_d, b))
    if abs(q) < 1e-16:
        return (-b / (2.0 * a),)
    return (q / a, c / q)


def search_level(level: int, include_all_gram_solutions: bool = False) -> Dict[str, object]:
    """Search solutions whose angles are integral multiples of pi/level.

    The symmetry inequalities are those in equation (6) of the paper.  The
    floating-point stage only proposes candidates; every returned tuple is
    certified exactly in a cyclotomic quotient and then checked geometrically.
    """

    if level < 3:
        raise ValueError("level must be at least 3")
    started = perf_counter()
    cosine = [cos(pi * m / level) for m in range(level + 1)]
    # The symmetry inequalities in equation (6) act on line angles.  Search
    # those internally, then convert theta to the interior dihedral angle
    # alpha=pi-theta before exact and geometric certification.
    proposed_line_angles: Set[Tuple[int, int, int, int, int, int]] = set()
    loop_count = 0

    for m12 in range(1, level):
        c12 = cosine[m12]
        for m34 in range(1, m12 + 1):
            c34 = cosine[m34]
            for m13 in range(1, level):
                c13 = cosine[m13]
                for m24 in range(1, m13 + 1):
                    if m13 + m24 > m12 + m34:
                        continue
                    c24 = cosine[m24]
                    for m14 in range(1, level):
                        loop_count += 1
                        c14 = cosine[m14]
                        A, B, C = quadratic_coefficients_for_c23(
                            (c12, c34, c13, c24, c14)
                        )
                        for root in _quadratic_roots(A, B, C):
                            if not isfinite(root) or root <= -1.0000000001 or root >= 1.0000000001:
                                continue
                            root = min(1.0, max(-1.0, root))
                            center = int(round(level * acos(root) / pi))
                            # Checking a small neighborhood makes the proposal
                            # stage insensitive to ordinary rounding error.
                            for m23 in range(max(1, center - 2), min(level - 1, center + 2) + 1):
                                if m14 + m23 > m13 + m24:
                                    continue
                                values = (m12, m34, m13, m24, m14, m23)
                                c = (c12, c34, c13, c24, c14, cosine[m23])
                                if abs(gram_determinant_from_cosines(c)) < 1e-9:
                                    proposed_line_angles.add(values)

    proposed = {
        tuple(level - m for m in line_values)
        for line_values in proposed_line_angles
    }
    exact = sorted(v for v in proposed if exact_torsion_zero(level, v))
    tetrahedra = [v for v in exact if is_nondegenerate_tetrahedron(level, v)]
    result: Dict[str, object] = {
        "level": level,
        "five_angle_loops": loop_count,
        "floating_candidates": len(proposed_line_angles),
        "exact_gram_solution_count": len(exact),
        "tetrahedron_count": len(tetrahedra),
        "tetrahedra": tetrahedra,
        "elapsed_seconds": perf_counter() - started,
    }
    if include_all_gram_solutions:
        result["exact_gram_solutions"] = exact
    return result
