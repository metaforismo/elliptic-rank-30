from __future__ import annotations

import json
import unittest
from pathlib import Path

from research.degree24_polynomial_section_obstruction import canonical_payload
from research.degree47_degenerate_obstruction import (
    build_certificate,
    verify_double_root_case,
    verify_local_order_classification,
    verify_triple_root_case,
)

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree47_degenerate_obstruction.json"


class Degree47DegenerateObstructionTests(unittest.TestCase):
    def test_local_order_classification(self) -> None:
        verify_local_order_classification()

    def test_triple_root_case(self) -> None:
        verify_triple_root_case()

    def test_double_root_case(self) -> None:
        verify_double_root_case()

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
