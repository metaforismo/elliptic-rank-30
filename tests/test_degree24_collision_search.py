from __future__ import annotations

import unittest
from fractions import Fraction

from research.degree24_collision_search import (
    exact_cube_root,
    parameter_values,
    phi,
    run_search,
)


class Degree24CollisionSearchTests(unittest.TestCase):
    def test_reduced_parameter_partition(self) -> None:
        values = list(parameter_values(2, 2))
        self.assertEqual(
            values,
            [
                Fraction(-2),
                Fraction(-1),
                Fraction(1),
                Fraction(2),
                Fraction(-1, 2),
                Fraction(1, 2),
            ],
        )
        self.assertEqual(len(values), len(set(values)))

    def test_phi_exact_value(self) -> None:
        self.assertEqual(phi(3), Fraction(-423, 256))

    def test_exact_cube_root(self) -> None:
        self.assertEqual(exact_cube_root(Fraction(8, 27)), Fraction(2, 3))
        self.assertEqual(exact_cube_root(Fraction(-125, 64)), Fraction(-5, 4))
        self.assertIsNone(exact_cube_root(Fraction(2, 3)))

    def test_small_box_has_no_collision(self) -> None:
        result = run_search(2, 2)
        self.assertEqual(result["parameter_count"], 6)
        self.assertEqual(result["collision_count"], 0)
        self.assertEqual(result["distinct_image_values"], 6)


if __name__ == "__main__":
    unittest.main()
