"""Independent exact-audit tools for rational dihedral tetrahedra."""

from .algebra import (
    gram_laurent_polynomial,
    exact_torsion_zero,
    verify_parametric_families,
)
from .geometry import is_nondegenerate_tetrahedron
from .known import SPORADIC_TETRAHEDRA

__all__ = [
    "SPORADIC_TETRAHEDRA",
    "exact_torsion_zero",
    "gram_laurent_polynomial",
    "is_nondegenerate_tetrahedron",
    "verify_parametric_families",
]

