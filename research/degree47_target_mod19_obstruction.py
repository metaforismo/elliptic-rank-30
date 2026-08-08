#!/usr/bin/env python3
"""Pure exact certificate excluding the generic degree-(4,7) target system.

The Wronskian reduction leaves three equations in rational b,d, with
c=2/11-3b.  After t=11b, z=121d, k=2-3t, exact combinations produce two
linear equations in z.  Their determinant is a primitive quintic P(t).
P has no root modulo 19 and its leading coefficient is nonzero modulo 19,
so P has no rational root.  The exceptional k=0 branch forces d=0 and is not
degree four.

Only Python standard-library exact arithmetic is used.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable

from research.degree24_polynomial_section_obstruction import (
    LaurentPolynomial,
    _variables,
    canonical_payload,
)


P_COEFFICIENTS_LOW = [
    -79764168,
    357341556,
    -533872086,
    266432683,
    -758384,
    154704,
]

P_MOD_19_VALUES = [17, 10, 8, 11, 1, 16, 15, 9, 2, 8, 11, 1, 10, 15, 12, 14, 11, 6, 13]


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return result


def polynomial_subtract(left: list[int], right: list[int]) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] -= right[index]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_value_mod(coefficients_low: Iterable[int], value: int, modulus: int) -> int:
    result = 0
    for coefficient in reversed(list(coefficients_low)):
        result = (result * value + coefficient) % modulus
    return result


def coefficient_gcd(coefficients: Iterable[int]) -> int:
    result = 0
    for coefficient in coefficients:
        result = math.gcd(result, abs(coefficient))
    return result


def verify_symbolic_reduction() -> None:
    t, z = _variables(2)
    k = 2 - 3 * t
    b = Fraction(1, 11) * t
    c = Fraction(1, 11) * k
    d = Fraction(1, 121) * z

    e8 = -9 * b**2 * c - 3 * b * c**2 + 36 * c * d + 32 * d**2
    e9 = (
        9 * b**4
        + 3 * b**3 * c
        - 72 * b**2 * d
        - 72 * b * c**2
        - 28 * b * c * d
        - 22 * c**3
        + 144 * d**2
    )
    e10 = (
        9 * b**3 * c
        + 3 * b**2 * c**2
        - 36 * b * c * d
        - 54 * c**3
        - 144 * c**2 * d
        + 256 * b**2 * d**2
    )

    h8 = 16 * z**2 + 198 * k * z - 33 * t * k
    if 11**4 * e8 != 2 * h8:
        raise AssertionError("scaled E8 identity failed")

    f9 = 11**4 * e9
    l9 = 12 * t**2 + 5290 * t - 3564
    c9 = -588 * t**3 - 4455 * t**2 + 6138 * t - 1936
    g9 = l9 * z + c9
    if f9 - 9 * h8 != g9:
        raise AssertionError("linear E9 combination failed")

    f10 = 11**6 * e10
    l10 = -264 * t**2 + 3993 * t - 2904
    c10 = 44 * t**3 - 53845 * t**2 + 71874 * t - 23958
    g10 = l10 * z + c10
    if f10 - 16 * t**2 * h8 != 12 * k * g10:
        raise AssertionError("linear E10 combination failed")

    determinant = c9 * l10 - c10 * l9
    expected = LaurentPolynomial(
        2,
        {
            (degree, 0): Fraction(coefficient)
            for degree, coefficient in enumerate(P_COEFFICIENTS_LOW)
        },
    )
    if determinant != expected:
        raise AssertionError("determinant quintic identity failed")


def verify_mod_19_obstruction() -> None:
    if coefficient_gcd(P_COEFFICIENTS_LOW) != 1:
        raise AssertionError("quintic is not primitive")
    if P_COEFFICIENTS_LOW[-1] % 19 == 0:
        raise AssertionError("leading coefficient vanishes modulo 19")
    values = [
        polynomial_value_mod(P_COEFFICIENTS_LOW, value, 19)
        for value in range(19)
    ]
    if values != P_MOD_19_VALUES:
        raise AssertionError(f"unexpected mod-19 values: {values}")
    if any(value == 0 for value in values):
        raise AssertionError("quintic has a root modulo 19")


def verify_determinant_by_integer_convolution() -> None:
    # Low-to-high coefficients of C9,L10,C10,L9.
    c9 = [-1936, 6138, -4455, -588]
    l10 = [-2904, 3993, -264]
    c10 = [-23958, 71874, -53845, 44]
    l9 = [-3564, 5290, 12]
    determinant = polynomial_subtract(
        polynomial_multiply(c9, l10),
        polynomial_multiply(c10, l9),
    )
    if determinant != P_COEFFICIENTS_LOW:
        raise AssertionError("integer-convolution determinant failed")


def build_certificate() -> dict:
    verify_symbolic_reduction()
    verify_determinant_by_integer_convolution()
    verify_mod_19_obstruction()
    return {
        "certificate_id": "degree47-target-mod19-obstruction-v1",
        "claim_status": "proved",
        "exceptional_branch": {
            "condition": "k=2-3t=0",
            "consequence": "the first reduced equation gives z=0, hence d=0",
            "result": "not a degree-four L",
        },
        "linear_equations": {
            "G9": "(12t^2+5290t-3564)z-588t^3-4455t^2+6138t-1936=0",
            "G10": "(-264t^2+3993t-2904)z+44t^3-53845t^2+71874t-23958=0",
        },
        "modular_obstruction": {
            "leading_coefficient_mod_19": P_COEFFICIENTS_LOW[-1] % 19,
            "modulus": 19,
            "values_at_0_through_18": P_MOD_19_VALUES,
            "zero_count": 0,
        },
        "normalization": {
            "c": "2/11-3b",
            "k": "2-3t",
            "t": "11b",
            "z": "121d",
        },
        "quintic": {
            "coefficients_low_to_high": P_COEFFICIENTS_LOW,
            "content": 1,
            "polynomial": "154704t^5-758384t^4+266432683t^3-533872086t^2+357341556t-79764168",
        },
        "rational_root_argument": "For a primitive integer polynomial, a rational root has denominator dividing the leading coefficient. Since 19 does not divide that coefficient, reduction gives a root in F_19, but the exhaustive value table has none.",
        "result": "no_rational_solution",
        "schema_version": 1,
        "scope": "generic squarefree/coprime degree-(4,7) target branch; repeated-root values are excluded separately",
        "theorem": "The generic degree-(4,7) polynomial target system has no rational solution.",
        "verification": {
            "determinant": "two independent exact checks: sparse Laurent-polynomial identity and integer convolution",
            "implementation": "Python standard library only",
            "modular_check": "complete enumeration of F_19",
            "scaled_equations": "exact sparse Laurent-polynomial arithmetic over Q",
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
    print("VERIFIED generic degree-(4,7) target obstruction modulo 19")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
