#!/usr/bin/env python3
"""Exact certificate for a restrictive polynomial-section ansatz.

The script proves that there are no rational coefficients with deg(L)=2 and
Q of degree at most 4 satisfying

    Q(v)^2 - v^2 L(v)^3 = v^3 - S v^2 + 3 v + 1.

Only exact rational and integer arithmetic from the Python standard library is
used.  This is an intermediate obstruction, not a rank-30 certificate.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

Exponent = Tuple[int, ...]


class LaurentPolynomial:
    """A tiny exact sparse Laurent-polynomial ring over Q."""

    def __init__(
        self,
        nvars: int,
        terms: Mapping[Exponent, Fraction | int] | None = None,
    ) -> None:
        self.nvars = nvars
        normalized: Dict[Exponent, Fraction] = {}
        for exponent, coefficient in (terms or {}).items():
            if len(exponent) != nvars:
                raise ValueError("wrong exponent length")
            value = Fraction(coefficient)
            if value:
                normalized[tuple(exponent)] = normalized.get(tuple(exponent), Fraction(0)) + value
        self.terms = {exponent: coefficient for exponent, coefficient in normalized.items() if coefficient}

    @classmethod
    def constant(cls, value: Fraction | int, nvars: int) -> "LaurentPolynomial":
        value = Fraction(value)
        if not value:
            return cls(nvars)
        return cls(nvars, {(0,) * nvars: value})

    @classmethod
    def variable(cls, index: int, nvars: int) -> "LaurentPolynomial":
        if not 0 <= index < nvars:
            raise IndexError(index)
        exponent = [0] * nvars
        exponent[index] = 1
        return cls(nvars, {tuple(exponent): Fraction(1)})

    def _coerce(self, other: object) -> "LaurentPolynomial":
        if isinstance(other, LaurentPolynomial):
            if other.nvars != self.nvars:
                raise ValueError("incompatible Laurent-polynomial rings")
            return other
        if isinstance(other, (int, Fraction)):
            return LaurentPolynomial.constant(other, self.nvars)
        return NotImplemented  # type: ignore[return-value]

    def __add__(self, other: object) -> "LaurentPolynomial":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        terms = dict(self.terms)
        for exponent, coefficient in rhs.terms.items():
            terms[exponent] = terms.get(exponent, Fraction(0)) + coefficient
        return LaurentPolynomial(self.nvars, terms)

    def __radd__(self, other: object) -> "LaurentPolynomial":
        return self + other

    def __neg__(self) -> "LaurentPolynomial":
        return LaurentPolynomial(self.nvars, {exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "LaurentPolynomial":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object) -> "LaurentPolynomial":
        lhs = self._coerce(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs - self

    def __mul__(self, other: object) -> "LaurentPolynomial":
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return NotImplemented
        terms: Dict[Exponent, Fraction] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in rhs.terms.items():
                exponent = tuple(a + b for a, b in zip(left_exponent, right_exponent))
                terms[exponent] = terms.get(exponent, Fraction(0)) + left_coefficient * right_coefficient
        return LaurentPolynomial(self.nvars, terms)

    def __rmul__(self, other: object) -> "LaurentPolynomial":
        return self * other

    def __pow__(self, exponent: int) -> "LaurentPolynomial":
        if not isinstance(exponent, int):
            raise TypeError("exponent must be an integer")
        if exponent < 0:
            if len(self.terms) != 1:
                raise ValueError("only monomials may have negative powers")
            monomial, coefficient = next(iter(self.terms.items()))
            return LaurentPolynomial(
                self.nvars,
                {tuple(exponent * value for value in monomial): coefficient**exponent},
            )
        result = LaurentPolynomial.constant(1, self.nvars)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power >>= 1
        return result

    def coefficient(self, variable_index: int, degree: int) -> "LaurentPolynomial":
        terms: Dict[Exponent, Fraction] = {}
        for exponent, coefficient in self.terms.items():
            if exponent[variable_index] != degree:
                continue
            reduced = list(exponent)
            reduced[variable_index] = 0
            key = tuple(reduced)
            terms[key] = terms.get(key, Fraction(0)) + coefficient
        return LaurentPolynomial(self.nvars, terms)

    def is_zero(self) -> bool:
        return not self.terms

    def __eq__(self, other: object) -> bool:
        rhs = self._coerce(other)
        if rhs is NotImplemented:
            return False
        return self.terms == rhs.terms


def _variables(nvars: int) -> Iterable[LaurentPolynomial]:
    return (LaurentPolynomial.variable(index, nvars) for index in range(nvars))


def verify_coefficient_extraction() -> None:
    """Recompute all nine coefficient equations formally and exactly."""

    a, b, c, d, e, f, g, h, s, v = _variables(10)
    linear = a * v**2 + b * v + c
    quartic = d * v**4 + e * v**3 + f * v**2 + g * v + h
    residual = quartic**2 - v**2 * linear**3 - (v**3 - s * v**2 + 3 * v + 1)

    expected = {
        0: h**2 - 1,
        1: 2 * g * h - 3,
        2: g**2 + 2 * f * h - c**3 + s,
        3: 2 * f * g + 2 * e * h - 3 * b * c**2 - 1,
        4: f**2 + 2 * e * g + 2 * d * h - 3 * a * c**2 - 3 * b**2 * c,
        5: 2 * d * g + 2 * e * f - 6 * a * b * c - b**3,
        6: e**2 + 2 * d * f - 3 * a**2 * c - 3 * a * b**2,
        7: 2 * d * e - 3 * a**2 * b,
        8: d**2 - a**3,
    }
    for degree, expression in expected.items():
        if residual.coefficient(9, degree) != expression:
            raise AssertionError(f"coefficient extraction failed in degree {degree}")


def verify_reduction_identities() -> None:
    """Verify the exact substitutions and the two remaining equations."""

    u, x, v = _variables(3)
    inverse_x = x**-1
    a = u**2
    b = u * x
    c = Fraction(1, 12) * (x**3 + 24) * inverse_x
    d = u**3
    e = Fraction(3, 2) * u**2 * x
    f = Fraction(1, 2) * u * (x**3 + 6) * inverse_x
    g = Fraction(3, 2)
    h = 1
    s = c**3 - Fraction(9, 4) - 2 * f

    linear = a * v**2 + b * v + c
    quartic = d * v**4 + e * v**3 + f * v**2 + g * v + h
    residual = quartic**2 - v**2 * linear**3 - (v**3 - s * v**2 + 3 * v + 1)

    for degree in (0, 1, 2, 5, 6, 7, 8):
        if not residual.coefficient(2, degree).is_zero():
            raise AssertionError(f"substitution failed in degree {degree}")

    square = (x**3 - 12) ** 2
    expected_degree_3 = 3 * u**2 * x - Fraction(1, 48) * u * square * x**-1 - 1
    expected_degree_4 = u**2 * (2 * u - Fraction(1, 48) * square * x**-2)
    if residual.coefficient(2, 3) != expected_degree_3:
        raise AssertionError("degree-3 residual identity failed")
    if residual.coefficient(2, 4) != expected_degree_4:
        raise AssertionError("degree-4 residual identity failed")

    first = 96 * u * x**2 - square
    second = 144 * u**2 * x**2 - u * square - 48 * x
    if 48 * x**2 * expected_degree_4 != u**2 * first:
        raise AssertionError("scaled degree-4 equation failed")
    if 48 * x * expected_degree_3 != second:
        raise AssertionError("scaled degree-3 equation failed")
    if second - u * first != 48 * x * (u**2 * x - 1):
        raise AssertionError("elimination identity failed")


def cubic_mod_7_values() -> list[int]:
    return [(24 * value**3 - 12 * value**2 - 6 * value - 1) % 7 for value in range(7)]


def valuation(value: Fraction, prime: int) -> int:
    if not value:
        raise ValueError("valuation of zero is infinite")
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def verify_final_obstruction() -> None:
    """Check the final factorization and the two rational obstructions."""

    # Coefficients are stored low degree first.
    expanded = [1, 0, -24, -96, 144]
    left = [-1, 6]
    right = [-1, -6, -12, 24]
    product = [0] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            product[i + j] += left_coefficient * right_coefficient
    if product != expanded:
        raise AssertionError("final quartic factorization failed")

    residues = cubic_mod_7_values()
    if residues != [6, 5, 5, 3, 3, 2, 4] or 0 in residues:
        raise AssertionError("the cubic unexpectedly has a root modulo 7")

    candidate = Fraction(1, 6)
    if valuation(candidate, 2) != -1:
        raise AssertionError("incorrect 2-adic valuation")
    if valuation(candidate, 2) % 3 == 0:
        raise AssertionError("1/6 was incorrectly accepted as a rational cube")


def build_certificate() -> dict:
    verify_coefficient_extraction()
    verify_reduction_identities()
    verify_final_obstruction()
    return {
        "certificate_id": "degree24-polynomial-section-obstruction-v1",
        "claim_status": "proved",
        "coefficient_equations": [
            "h^2 = 1",
            "2 g h = 3",
            "g^2 + 2 f h - c^3 = -S",
            "2 f g + 2 e h - 3 b c^2 = 1",
            "f^2 + 2 e g + 2 d h = 3 a c^2 + 3 b^2 c",
            "2 d g + 2 e f = 6 a b c + b^3",
            "e^2 + 2 d f = 3 a^2 c + 3 a b^2",
            "2 d e = 3 a^2 b",
            "d^2 = a^3",
        ],
        "cube_obstruction": {
            "candidate": "z = 1/6",
            "conclusion": "not a cube in Q",
            "prime": 2,
            "valuation": -1,
            "valuation_mod_3": 2,
        },
        "elimination": {
            "D": "(x^3 - 12)^2",
            "deduction": "u^2 x = 1",
            "factorization": "(6 z - 1)(24 z^3 - 12 z^2 - 6 z - 1)",
            "final_equation": "144 z^4 - 96 z^3 - 24 z^2 + 1 = 0",
            "linear_combination": "second - u*first = 48 x (u^2 x - 1)",
            "nonzero_conditions": ["u != 0", "x != 0"],
            "remaining_equations": [
                "D = 96 u x^2",
                "144 u^2 x^2 - u D - 48 x = 0",
            ],
            "z_substitution": "z = u^3",
        },
        "identity": "Q(v)^2 - v^2 L(v)^3 = v^3 - S v^2 + 3 v + 1",
        "modular_irreducibility": {
            "conclusion": "the cubic has no rational root",
            "cubic": "24 z^3 - 12 z^2 - 6 z - 1",
            "has_root_mod_7": False,
            "modulus": 7,
            "values": [6, 5, 5, 3, 3, 2, 4],
        },
        "normalization": {
            "h": "1 after replacing Q by -Q if necessary",
            "g": "3/2",
            "leading_parameter": "u = d/a, hence a = u^2 and d = u^3",
        },
        "result": "no_rational_solution",
        "schema_version": 1,
        "scope": {
            "L": "a v^2 + b v + c with a != 0",
            "Q": "d v^4 + e v^3 + f v^2 + g v + h",
            "field": "Q",
            "limitation": "This rules out only this degree-(2,4) polynomial-section ansatz; it is not a rank-30 obstruction.",
        },
        "substitutions": {
            "b": "u x",
            "c": "(x^3 + 24)/(12 x)",
            "e": "(3/2) u^2 x",
            "f": "u (x^3 + 6)/(2 x)",
        },
        "theorem": "No rational coefficients with deg(L)=2 and deg(Q)<=4 satisfy the stated identity for any S in Q.",
        "verification": {
            "coefficient_extraction": "exact sparse Laurent-polynomial arithmetic",
            "factorization": "exact integer convolution",
            "implementation": "Python standard library only",
            "modular_root_check": "complete enumeration in F_7",
            "reduction_identities": "exact sparse Laurent-polynomial arithmetic",
        },
    }


def canonical_payload(certificate: dict) -> str:
    return json.dumps(certificate, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the canonical certificate JSON")
    parser.add_argument("--check", type=Path, help="compare against an existing certificate JSON")
    args = parser.parse_args()

    payload = canonical_payload(build_certificate())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check:
        existing = args.check.read_text(encoding="utf-8")
        if existing != payload:
            raise SystemExit(f"certificate mismatch: {args.check}")
    print("VERIFIED degree-(2,4) polynomial-section obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
