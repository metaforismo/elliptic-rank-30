from __future__ import annotations

import json
import unittest
from pathlib import Path

from research.degree24_polynomial_section_obstruction import (
    build_certificate,
    canonical_payload,
    cubic_mod_7_values,
    verify_coefficient_extraction,
    verify_final_obstruction,
    verify_reduction_identities,
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree24_polynomial_section_obstruction.json"


class Degree24PolynomialSectionObstructionTests(unittest.TestCase):
    def test_general_coefficient_extraction(self) -> None:
        verify_coefficient_extraction()

    def test_reduction_identities(self) -> None:
        verify_reduction_identities()

    def test_final_obstruction(self) -> None:
        verify_final_obstruction()
        self.assertEqual(cubic_mod_7_values(), [6, 5, 5, 3, 3, 2, 4])

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
