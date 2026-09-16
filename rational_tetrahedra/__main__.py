"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
from time import perf_counter

from .algebra import exact_torsion_zero, gram_laurent_polynomial, integral_gram_terms, verify_parametric_families
from .geometry import is_nondegenerate_tetrahedron
from .known import SPORADIC_TETRAHEDRA
from .search import search_level


def audit() -> dict:
    started = perf_counter()
    poly = gram_laurent_polynomial()
    scale, terms = integral_gram_terms()
    exact_failures = []
    geometry_failures = []
    for level, values in SPORADIC_TETRAHEDRA:
        if not exact_torsion_zero(level, values):
            exact_failures.append((level, values))
        if not is_nondegenerate_tetrahedron(level, values):
            geometry_failures.append((level, values))
    return {
        "laurent_monomials": len(poly),
        "integral_scale": scale,
        "integral_constant_coefficient": dict(terms).get((0, 0, 0, 0, 0, 0)),
        "parametric_families": verify_parametric_families(),
        "sporadic_count": len(SPORADIC_TETRAHEDRA),
        "sporadic_exact_failures": exact_failures,
        "sporadic_geometry_failures": geometry_failures,
        "elapsed_seconds": perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("audit", help="audit the polynomial, families, and Table 3")
    search_parser = sub.add_parser("search", help="run the bounded low-order search")
    search_parser.add_argument("level", type=int)
    search_parser.add_argument(
        "--include-all-gram-solutions",
        action="store_true",
        help="also print every certified Gram solution, including non-tetrahedra",
    )
    args = parser.parse_args()

    result = (
        audit()
        if args.command == "audit"
        else search_level(args.level, args.include_all_gram_solutions)
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
