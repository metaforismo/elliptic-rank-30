#!/usr/bin/env python3
"""Exact classification of a degree-(2,4) polynomial decomposition.

For p != 0, classify all rational solutions of

    Q(v)^2 - v^2 L(v)^3 = v^3 - S v^2 + p v + 1

with deg(L)=2 and deg(Q)<=4.  The classification is a one-parameter rational
family subject to a cube-class condition.  This is an intermediate construction
result, not a rank-30 certificate.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict

from research.degree24_polynomial_section_obstruction import (
    LaurentPolynomial,
    _variables,
    canonical_payload,
)


def verify_general_reduction() -> None:
    """Recompute the reduced coefficient scheme with p left symbolic."""

    u, x, p, v = _variables(4)
    c = Fraction(1, 12) * (x**3 + 8 * p) * x**-1
    a = u**2
    b = u * x
    d = u**3
    e = Fraction(3, 2) * u**2 * x
    f = Fraction(1, 2) * u * (x**3 + 2 * p) * x**-1
    g = Fraction(1, 2) * p
    h = 1
    s = c**3 - g**2 - 2 * f

    linear = a * v**2 + b * v + c
    quartic = d * v**4 + e * v**3 + f * v**2 + g * v + h
    residual = quartic**2 - v**2 * linear**3 - (v**3 - s * v**2 + p * v + 1)

    for degree in (0, 1, 2, 5, 6, 7, 8):
        if not residual.coefficient(3, degree).is_zero():
            raise AssertionError(f"general reduction failed in degree {degree}")

    square = (x**3 - 4 * p) ** 2
    expected_degree_3 = 3 * u**2 * x - Fraction(1, 48) * u * square * x**-1 - 1
    expected_degree_4 = u**2 * (2 * u - Fraction(1, 48) * square * x**-2)
    if residual.coefficient(3, 3) != expected_degree_3:
        raise AssertionError("general degree-3 identity failed")
    if residual.coefficient(3, 4) != expected_degree_4:
        raise AssertionError("general degree-4 identity failed")

    first = 96 * u * x**2 - square
    second = 144 * u**2 * x**2 - u * square - 48 * x
    if second - u * first != 48 * x * (u**2 * x - 1):
        raise AssertionError("general elimination identity failed")


def verify_parameterization_formally() -> None:
    """Substitute the claimed t-family and verify the full identity exactly."""

    t, v = _variables(2)
    p = Fraction(1, 186624) * t**12 - Fraction(1, 6) * t**3
    a = 36 * t**-4
    b = Fraction(1, 6) * t**2
    c = Fraction(1, 5184) * t**8 - 4 * t**-1
    d = 216 * t**-6
    e = Fraction(3, 2)
    f = Fraction(1, 288) * t**6 - 36 * t**-3
    g = Fraction(1, 2) * p
    h = 1
    s = c**3 - g**2 - 2 * f

    linear = a * v**2 + b * v + c
    quartic = d * v**4 + e * v**3 + f * v**2 + g * v + h
    residual = quartic**2 - v**2 * linear**3 - (v**3 - s * v**2 + p * v + 1)
    if not residual.is_zero():
        raise AssertionError("the t-parameterization does not satisfy the identity")

    u = 6 * t**-2
    x = Fraction(1, 36) * t**4
    square = (x**3 - 4 * p) ** 2
    if u**2 * x != 1:
        raise AssertionError("u^2 x parameter identity failed")
    if square != 96 * u * x**2:
        raise AssertionError("square parameter identity failed")


def coefficients(t: Fraction | int) -> Dict[str, Fraction]:
    t = Fraction(t)
    if not t:
        raise ValueError("t must be nonzero")
    p = t**12 / 186624 - t**3 / 6
    a = 36 / t**4
    b = t**2 / 6
    c = t**8 / 5184 - 4 / t
    d = 216 / t**6
    e = Fraction(3, 2)
    f = t**6 / 288 - 36 / t**3
    g = p / 2
    h = Fraction(1)
    s = c**3 - g**2 - 2 * f
    return {
        "p": p,
        "S": s,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "e": e,
        "f": f,
        "g": g,
        "h": h,
    }


def _poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    result = [Fraction(0)] * size
    for index in range(size):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] += right[index]
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def _poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def _poly_pow(base: list[Fraction], exponent: int) -> list[Fraction]:
    result = [Fraction(1)]
    power = list(base)
    while exponent:
        if exponent & 1:
            result = _poly_mul(result, power)
        power = _poly_mul(power, power)
        exponent >>= 1
    return result


def verify_numeric_identity(t: Fraction | int) -> None:
    values = coefficients(t)
    linear = [values["c"], values["b"], values["a"]]
    quartic = [values["h"], values["g"], values["f"], values["e"], values["d"]]
    left = _poly_add(
        _poly_pow(quartic, 2),
        [-value for value in ([Fraction(0), Fraction(0)] + _poly_pow(linear, 3))],
    )
    right = [Fraction(1), values["p"], -values["S"], Fraction(1)]
    if left != right:
        raise AssertionError(f"numeric identity failed at t={t}")


def verify_target_p3_obstruction() -> None:
    """Use q=t^3/36 to prove p=3 is outside the rational parameter image."""

    # 3 q^4 - 2 q - 1 = (q-1)(3 q^3 + 3 q^2 + 3 q + 1).
    left = [Fraction(-1), Fraction(-2), Fraction(0), Fraction(0), Fraction(3)]
    factor_one = [Fraction(-1), Fraction(1)]
    factor_two = [Fraction(1), Fraction(3), Fraction(3), Fraction(3)]
    if _poly_mul(factor_one, factor_two) != left:
        raise AssertionError("p=3 factorization failed")

    candidates = [Fraction(1), Fraction(-1), Fraction(1, 3), Fraction(-1, 3)]
    cubic_values = [3 * q**3 + 3 * q**2 + 3 * q + 1 for q in candidates]
    if any(value == 0 for value in cubic_values):
        raise AssertionError("the residual cubic has a rational root")

    # The sole rational root of the quartic is q=1.  It would require t^3=36,
    # impossible because v_2(36)=2 is not divisible by 3.
    numerator = 36
    exponent_two = 0
    while numerator % 2 == 0:
        numerator //= 2
        exponent_two += 1
    if exponent_two != 2 or exponent_two % 3 == 0:
        raise AssertionError("cube-class obstruction for t^3=36 failed")


def _fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def build_certificate() -> dict:
    verify_general_reduction()
    verify_parameterization_formally()
    for sample in (Fraction(1), Fraction(2), Fraction(3), Fraction(-3), Fraction(5, 2)):
        verify_numeric_identity(sample)
    verify_target_p3_obstruction()

    sample = coefficients(3)
    return {
        "certificate_id": "degree24-polynomial-section-classification-v1",
        "claim_status": "proved",
        "classification": {
            "assumption": "p != 0",
            "if_and_only_if": "a rational solution exists iff p=t^12/186624-t^3/6 for some t in Q^*",
            "parameter": "t in Q^*",
            "parameter_map": "p(t)=t^12/186624-t^3/6",
        },
        "coefficient_formulas": {
            "L": {
                "a": "36/t^4",
                "b": "t^2/6",
                "c": "t^8/5184-4/t",
            },
            "Q": {
                "d": "216/t^6",
                "e": "3/2",
                "f": "t^6/288-36/t^3",
                "g": "p(t)/2",
                "h": "1",
            },
            "S": "c^3-g^2-2f",
        },
        "identity": "Q(v)^2-v^2 L(v)^3=v^3-S v^2+p v+1",
        "inverse_construction_view": {
            "cube_coset_variable": "q=t^3/36",
            "quartic_map": "p=9q^4-6q",
            "research_use": "search high-multiplicity fibers and combine independent forcing mechanisms",
        },
        "p_equals_3": {
            "cube_class_obstruction": "q=1 would require t^3=36, but v_2(36)=2 is not divisible by 3",
            "factorization": "3q^4-2q-1=(q-1)(3q^3+3q^2+3q+1)",
            "residual_cubic_rational_roots": [],
            "result": "no rational degree-(2,4) decomposition",
        },
        "sample_t_3": {key: _fraction_string(value) for key, value in sample.items()},
        "schema_version": 1,
        "scope_limitation": "This classifies only the stated polynomial identity and does not itself construct a rank-30 elliptic curve.",
        "theorem": "For p nonzero, the complete rational degree-(2,4) coefficient scheme is the displayed one-parameter family.",
        "verification": {
            "forward_reduction": "exact symbolic Laurent-polynomial identities",
            "reverse_substitution": "exact symbolic Laurent-polynomial identity",
            "sample_checks": ["1", "2", "3", "-3", "5/2"],
            "target_obstruction": "exact factorization, rational-root theorem, and 2-adic cube valuation",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()

    payload = canonical_payload(build_certificate())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check and args.check.read_text(encoding="utf-8") != payload:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print("VERIFIED complete degree-(2,4) classification")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
