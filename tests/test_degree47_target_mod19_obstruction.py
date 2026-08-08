from __future__ import annotations

import json
import unittest
from pathlib import Path

from research.degree24_polynomial_section_obstruction import canonical_payload
from research.degree47_target_mod19_obstruction import (
    build_certificate,
    verify_determinant_by_integer_convolution,
    verify_mod_19_obstruction,
    verify_symbolic_reduction,
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree47_target_mod19_obstruction.json"


class Degree47TargetMod19ObstructionTests(unittest.TestCase):
    def test_symbolic_reduction(self) -> None:
        verify_symbolic_reduction()

    def test_integer_determinant(self) -> None:
        verify_determinant_by_integer_convolution()

    def test_mod_19_obstruction(self) -> None:
        verify_mod_19_obstruction()

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
