from __future__ import annotations

import json
import unittest
from fractions import Fraction
from pathlib import Path

from research.degree24_polynomial_section_classification import (
    build_certificate,
    coefficients,
    verify_general_reduction,
    verify_numeric_identity,
    verify_parameterization_formally,
    verify_target_p3_obstruction,
)
from research.degree24_polynomial_section_obstruction import canonical_payload

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree24_polynomial_section_classification.json"


class Degree24PolynomialSectionClassificationTests(unittest.TestCase):
    def test_forward_reduction(self) -> None:
        verify_general_reduction()

    def test_reverse_parameterization(self) -> None:
        verify_parameterization_formally()

    def test_exact_samples(self) -> None:
        for parameter in (Fraction(1), Fraction(2), Fraction(3), Fraction(-3), Fraction(5, 2)):
            verify_numeric_identity(parameter)

    def test_t3_sample(self) -> None:
        self.assertEqual(
            coefficients(3),
            {
                "p": Fraction(-423, 256),
                "S": Fraction(-665, 216),
                "a": Fraction(4, 9),
                "b": Fraction(3, 2),
                "c": Fraction(-13, 192),
                "d": Fraction(8, 27),
                "e": Fraction(3, 2),
                "f": Fraction(115, 96),
                "g": Fraction(-423, 512),
                "h": Fraction(1),
            },
        )

    def test_target_p3_obstruction(self) -> None:
        verify_target_p3_obstruction()

    def test_committed_certificate_matches_recomputation(self) -> None:
        expected = build_certificate()
        observed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(CERTIFICATE.read_text(encoding="utf-8"), canonical_payload(expected))


if __name__ == "__main__":
    unittest.main()
