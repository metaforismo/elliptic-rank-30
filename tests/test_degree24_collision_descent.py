from __future__ import annotations

import json
import unittest
from pathlib import Path

from research.degree24_collision_descent import (
    build_certificate,
    verify_gcd_identities,
    verify_three_adic_residues,
    verify_two_adic_residues,
)
from research.degree24_polynomial_section_obstruction import canonical_payload

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree24_collision_descent.json"


class Degree24CollisionDescentTests(unittest.TestCase):
    def test_gcd_identities(self) -> None:
        verify_gcd_identities()

    def test_two_adic_residues(self) -> None:
        verify_two_adic_residues()

    def test_three_adic_residues(self) -> None:
        verify_three_adic_residues()

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
