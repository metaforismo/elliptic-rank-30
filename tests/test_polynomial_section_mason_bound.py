from __future__ import annotations

import json
import unittest
from pathlib import Path

from research.degree24_polynomial_section_obstruction import canonical_payload
from research.polynomial_section_mason_bound import (
    build_certificate,
    maximal_term_degree,
    mason_upper_radical_degree,
    verify_cubic_degeneracies,
    verify_mason_inequality,
    verify_passport,
    verify_wronskian_degree,
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "polynomial_section_mason_bound.json"


class PolynomialSectionMasonBoundTests(unittest.TestCase):
    def test_mason_inequality(self) -> None:
        verify_mason_inequality()
        self.assertEqual(maximal_term_degree(2, 0), 14)
        self.assertEqual(mason_upper_radical_degree(2, 0), 15)
        self.assertGreater(maximal_term_degree(3, 0), mason_upper_radical_degree(3, 0) - 1)

    def test_passport(self) -> None:
        verify_passport()

    def test_target_degeneracies(self) -> None:
        verify_cubic_degeneracies()

    def test_wronskian_degree(self) -> None:
        verify_wronskian_degree()

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
