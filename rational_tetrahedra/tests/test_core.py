from __future__ import annotations

import unittest

from rational_tetrahedra.algebra import (
    exact_torsion_zero,
    gram_laurent_polynomial,
    integral_gram_terms,
    verify_parametric_families,
)
from rational_tetrahedra.geometry import is_nondegenerate_tetrahedron
from rational_tetrahedra.known import SPORADIC_TETRAHEDRA
from rational_tetrahedra.search import search_level


class RationalTetrahedraAuditTests(unittest.TestCase):
    def test_gram_laurent_shape(self) -> None:
        poly = gram_laurent_polynomial()
        scale, terms = integral_gram_terms()
        self.assertEqual(len(poly), 105)
        self.assertEqual(scale, 16)
        self.assertEqual(dict(terms)[(0, 0, 0, 0, 0, 0)], -20)

    def test_parametric_families(self) -> None:
        self.assertEqual(
            verify_parametric_families(),
            {"family_1": True, "family_2": True},
        )

    def test_table_3(self) -> None:
        self.assertEqual(len(SPORADIC_TETRAHEDRA), 59)
        for level, values in SPORADIC_TETRAHEDRA:
            with self.subTest(level=level, values=values):
                self.assertTrue(exact_torsion_zero(level, values))
                self.assertTrue(is_nondegenerate_tetrahedron(level, values))

    def test_level_12_search(self) -> None:
        result = search_level(12)
        self.assertEqual(result["tetrahedron_count"], 7)
        found = {tuple(v) for v in result["tetrahedra"]}
        self.assertIn((3, 4, 3, 4, 6, 8), found)
        self.assertIn((3, 6, 4, 6, 4, 6), found)


if __name__ == "__main__":
    unittest.main()
