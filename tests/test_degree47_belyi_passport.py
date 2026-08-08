from __future__ import annotations

import unittest

from research.degree47_belyi_passport import (
    canonical_sigma_infinity,
    centralizer_of_sigma_infinity,
    compose,
    cycle_partition,
    identity,
    inverse,
    is_transitive,
    perfect_matchings,
)


class Degree47BelyiPassportTests(unittest.TestCase):
    def test_permutation_arithmetic(self) -> None:
        permutation = (1, 2, 0, 3)
        self.assertEqual(compose(permutation, inverse(permutation)), identity(4))
        self.assertEqual(cycle_partition(permutation), (3, 1))

    def test_fixed_infinity_and_centralizer(self) -> None:
        sigma_infinity = canonical_sigma_infinity()
        self.assertEqual(cycle_partition(sigma_infinity), (11, 1, 1, 1))
        centralizer = centralizer_of_sigma_infinity()
        self.assertEqual(len(centralizer), 66)
        self.assertTrue(
            all(
                compose(element, sigma_infinity) == compose(sigma_infinity, element)
                for element in centralizer
            )
        )

    def test_small_perfect_matching_count(self) -> None:
        matchings = list(perfect_matchings(tuple(range(6))))
        self.assertEqual(len(matchings), 15)
        self.assertEqual(len(set(matchings)), 15)
        self.assertTrue(all(cycle_partition(matching) == (2, 2, 2) for matching in matchings))

    def test_transitivity(self) -> None:
        cycle = (1, 2, 3, 0)
        self.assertTrue(is_transitive((cycle,)))
        self.assertFalse(is_transitive(((1, 0, 3, 2),)))


if __name__ == "__main__":
    unittest.main()
