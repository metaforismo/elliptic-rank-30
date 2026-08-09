from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research"
if str(RESEARCH) not in sys.path:
    sys.path.insert(0, str(RESEARCH))

from certify_three_trisection_plane_cubic import KNOWN_ROOTS, KNOWN_T, ORIGIN, packet_points_mod
from plane_cubic_finite_reduction import PlaneCubicGroup, row_rank


class PlaneCubicFiniteReductionTests(unittest.TestCase):
    def test_known_group_orders(self) -> None:
        expected = {43: 52, 53: 60, 61: 76, 67: 73, 101: 104, 191: 192}
        for prime, order in expected.items():
            group = PlaneCubicGroup.create(prime, KNOWN_T.numerator, KNOWN_T.denominator, ORIGIN)
            self.assertEqual(len(group.points), order)

    def test_known_packet_has_mod2_rank_ten(self) -> None:
        rows = [[] for _ in range(14)]
        torsion_witness = None
        for prime in [61, 67, 101, 179, 211, 263]:
            try:
                group = PlaneCubicGroup.create(prime, KNOWN_T.numerator, KNOWN_T.denominator, ORIGIN)
                points, labels, _ = packet_points_mod(group, KNOWN_ROOTS)
                vectors, dimension = group.quotient_vectors_points(points, 2)
            except (AssertionError, ValueError, ZeroDivisionError):
                continue
            self.assertEqual(len(labels), 14)
            if len(group.points) % 2:
                torsion_witness = prime
            if dimension:
                for index, vector in enumerate(vectors):
                    rows[index].extend(vector)
        self.assertEqual(row_rank(rows, 2), 10)
        self.assertEqual(torsion_witness, 67)

    def test_inverse_and_identity_on_nonflex_origin(self) -> None:
        group = PlaneCubicGroup.create(61, KNOWN_T.numerator, KNOWN_T.denominator, ORIGIN)
        for point in group.points:
            self.assertEqual(group.add(point, group.origin), point)
            self.assertEqual(group.add(point, group.neg(point)), group.origin)


if __name__ == "__main__":
    unittest.main()
