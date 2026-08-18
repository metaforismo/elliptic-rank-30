#!/usr/bin/env sage -python
"""Compatibility launcher for the open three-variable IV elimination.

Sage 10.9's multivariate libSingular polynomial class does not expose
``quo_rem``. The mathematical implementation is kept unchanged except for
replacing the attempted exact division by the determinant with an equivalent
fraction-field denominator test.
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "rank17_iv_three_variable_open_elimination_sage.py"

OLD = """        while True:\n            quotient, remainder = numerator.quo_rem(expected_determinant)\n            if remainder:\n                break\n            numerator = quotient\n            determinant_power += 1\n"""

NEW = """        while True:\n            candidate = fraction(numerator) / fraction(expected_determinant)\n            if candidate.denominator() != 1:\n                break\n            numerator = base(candidate.numerator())\n            determinant_power += 1\n"""


def main() -> None:
    source = IMPLEMENTATION.read_text(encoding="utf-8")
    if source.count(OLD) != 1:
        raise AssertionError("the expected Sage-10.9 compatibility site changed")
    patched = source.replace(OLD, NEW)
    namespace = {
        "__name__": "rank17_iv_three_variable_open_elimination_sage_v2_body",
        "__file__": str(IMPLEMENTATION),
    }
    exec(compile(patched, str(IMPLEMENTATION), "exec"), namespace)
    namespace["main"]()


if __name__ == "__main__":
    main()
