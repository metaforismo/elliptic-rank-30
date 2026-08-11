#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

PARAMS = {
    "e0": 1,
    "e1": 1,
    "p0": -3,
    "p1": 32,
    "p2": 66,
    "p3": 768,
    "q0": -12,
    "q1": -12,
    "q2": 0,
    "s": 0,
}
EXPECTED_C4 = [9, -192, 628, -384, 30, 64, -12, 0, 1]
EXPECTED_C6 = [
    -27,
    864,
    -7434,
    15488,
    -14841,
    2880,
    2316,
    -1152,
    99,
    96,
    -18,
    0,
    1,
]
EXPECTED_QUARTIC = [-3, 56, -6, 0, 1]


def trim(values):
    result = list(values)
    while result and result[-1] == 0:
        result.pop()
    return result


def add(left, right):
    size = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        ]
    )


def subtract(left, right):
    return add(left, [-value for value in right])


def multiply(left, right):
    if not left or not right:
        return []
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return trim(result)


def scale(values, scalar):
    return trim([scalar * value for value in values])


def power(values, exponent):
    result = [1]
    base = values
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def divide_exact(dividend, divisor):
    remainder = list(dividend)
    divisor = trim(divisor)
    quotient = [0] * max(0, len(remainder) - len(divisor) + 1)
    while len(trim(remainder)) >= len(divisor):
        remainder = trim(remainder)
        shift = len(remainder) - len(divisor)
        if remainder[-1] % divisor[-1]:
            raise AssertionError("nonintegral polynomial division")
        coefficient = remainder[-1] // divisor[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
    if trim(remainder):
        raise AssertionError(("nonzero remainder", remainder))
    return trim(quotient)


def gcd_over_q(left, right):
    left = [Fraction(value) for value in trim(left)]
    right = [Fraction(value) for value in trim(right)]
    while right:
        remainder = left[:]
        while len(remainder) >= len(right):
            coefficient = remainder[-1] / right[-1]
            shift = len(remainder) - len(right)
            for index, value in enumerate(right):
                remainder[index + shift] -= coefficient * value
            while remainder and remainder[-1] == 0:
                remainder.pop()
        left, right = right, remainder
    if not left:
        return []
    leading = left[-1]
    return [value / leading for value in left]


def derivative(values):
    return [index * values[index] for index in range(1, len(values))]


def evaluate(values, argument):
    result = 0
    for coefficient in reversed(values):
        result = result * argument + coefficient
    return result


def build_invariants():
    p0 = PARAMS["p0"]
    p1 = PARAMS["p1"]
    p2 = PARAMS["p2"]
    p3 = PARAMS["p3"]
    q0 = PARAMS["q0"]
    q1 = PARAMS["q1"]
    q2 = PARAMS["q2"]
    s = PARAMS["s"]
    e0 = PARAMS["e0"]
    e1 = PARAMS["e1"]

    c4 = [0] * 9
    c4[0] = p0 * p0
    c4[1] = 2 * p0 * p1
    c4[2] = 2 * p0 * p2 + p1 * p1
    c4[3] = 2 * p0 * p3 + 2 * p1 * p2
    c4[4] = (
        -3
        - 15 * p0 * p0
        - 20 * p0 * p1
        - 12 * p0 * p2
        - 6 * p0 * p3
        - 6 * p1 * p1
        - 6 * p1 * p2
        + 15 * q0 * q0
        - 10 * q0 * q1
        + 2 * q0 * q2
        + q1 * q1
        - s
    )
    c4[5] = (
        8
        + 24 * p0 * p0
        + 30 * p0 * p1
        + 16 * p0 * p2
        + 6 * p0 * p3
        + 8 * p1 * p1
        + 6 * p1 * p2
        - 24 * q0 * q0
        + 18 * q0 * q1
        - 4 * q0 * q2
        - 2 * q1 * q1
        + 3 * s
    )
    c4[6] = (
        -6
        - 10 * p0 * p0
        - 12 * p0 * p1
        - 6 * p0 * p2
        - 2 * p0 * p3
        - 3 * p1 * p1
        - 2 * p1 * p2
        + 10 * q0 * q0
        - 8 * q0 * q1
        + 2 * q0 * q2
        + q1 * q1
        - 3 * s
    )
    c4[7] = s
    c4[8] = 1

    reversed_c4 = list(reversed(c4))
    cubed = power(reversed_c4, 3)
    square_root = [0] * 12
    square_root[0] = 1
    for order in range(1, 12):
        numerator = cubed[order] - sum(
            square_root[index] * square_root[order - index]
            for index in range(1, order)
        )
        if numerator % 2:
            raise AssertionError(("odd square-root recurrence", order))
        square_root[order] = numerator // 2

    c6 = [0] * 13
    c6[12] = 1
    for index in range(1, 12):
        c6[index] = square_root[12 - index]
    c6[0] = e0 * p0**3

    equations = [
        c6[1] - 3 * e0 * p0 * p0 * p1,
        c6[2] - 3 * e0 * (p0 * p0 * p2 + p0 * p1 * p1),
        c6[3]
        - e0
        * (3 * p0 * p0 * p3 + 6 * p0 * p1 * p2 + p1**3),
        sum(c6) - e1 * q0**3,
        sum(index * c6[index] for index in range(13))
        - 3 * e1 * q0 * q0 * q1,
        sum(
            (index * (index - 1) // 2) * c6[index]
            for index in range(13)
        )
        - 3 * e1 * (q0 * q0 * q2 + q0 * q1 * q1),
    ]
    return c4, c6, equations


def compute_certificate():
    c4, c6, equations = build_invariants()
    if c4 != EXPECTED_C4 or c6 != EXPECTED_C6:
        raise AssertionError("unexpected invariant polynomials")
    if equations != [0] * 6:
        raise AssertionError(("Hermite--Pade equations failed", equations))

    numerator = subtract(power(c4, 3), power(c6, 2))
    fixed_factor = multiply(power([0, 1], 4), power([-1, 1], 3))
    residual = divide_exact(numerator, fixed_factor)
    expected_residual = scale(
        multiply([3, 1], EXPECTED_QUARTIC), 139968
    )
    if residual != expected_residual:
        raise AssertionError("unexpected discriminant factorization")
    if gcd_over_q(EXPECTED_QUARTIC, derivative(EXPECTED_QUARTIC)) != [
        Fraction(1)
    ]:
        raise AssertionError("quartic factor is not squarefree")
    if gcd_over_q(c4, multiply([3, 1], EXPECTED_QUARTIC)) != [
        Fraction(1)
    ]:
        raise AssertionError("residual discriminant meets c4")

    split0 = -3 * PARAMS["e0"] * PARAMS["p0"]
    split1 = -3 * PARAMS["e1"] * PARAMS["q0"]
    if math.isqrt(split0) ** 2 != split0:
        raise AssertionError("I4 tangent cone does not split")
    if math.isqrt(split1) ** 2 != split1:
        raise AssertionError("I3 tangent cone does not split")

    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_exact_split_semistable_seed_model",
        "truth_status": (
            "CERTIFIED exact algebraic model; Mordell-Weil rank-one "
            "section not yet constructed"
        ),
        "parameters": PARAMS,
        "normalized_invariants": {
            "c4_coefficients_ascending": c4,
            "c6_coefficients_ascending": c6,
        },
        "short_weierstrass_equation": (
            "y^2=x^3-3*c4(t)*x-2*c6(t)"
        ),
        "six_hermite_pade_equations": [str(value) for value in equations],
        "discriminant_numerator": {
            "definition": "c4(t)^3-c6(t)^2",
            "factorization": (
                "139968*t^4*(t-1)^3*(t+3)*"
                "(t^4-6*t^2+56*t-3)"
            ),
            "residual_quintic_coefficients_ascending": residual,
            "residual_squarefree": True,
            "coprime_to_c4": True,
        },
        "fiber_data": {
            "t=0": {
                "order": 4,
                "c4_value": evaluate(c4, 0),
                "type": "split I4",
                "tangent_square_coefficient": split0,
                "tangent_roots": [
                    -math.isqrt(split0),
                    math.isqrt(split0),
                ],
            },
            "t=1": {
                "order": 3,
                "c4_value": evaluate(c4, 1),
                "type": "split I3",
                "tangent_square_coefficient": split1,
                "tangent_roots": [
                    -math.isqrt(split1),
                    math.isqrt(split1),
                ],
            },
            "t=infinity": {
                "order": 12,
                "leading_c4": c4[-1],
                "type": (
                    "I12 over algebraic closure; splitness not asserted"
                ),
            },
            "additional_fibers": (
                "five simple I1 fibers over the roots of "
                "(t+3)*(t^4-6*t^2+56*t-3)"
            ),
        },
        "cross_prime_reconstruction": {
            "fixed_parameters": {
                "p0": -3,
                "q0": -12,
                "e0": 1,
                "e1": 1,
            },
            "variables": ["p1", "p2", "p3", "q1", "q2", "s"],
            "mod_5_seed": [2, 1, 3, 3, 0, 0],
            "mod_7_seed": [4, 3, 5, 2, 0, 0],
            "fixed_variable_jacobian_determinants": {"5": 2, "7": 6},
            "common_exact_reconstruction": [32, 66, 768, -12, 0, 0],
        },
        "limitations": [
            (
                "The model satisfies the required split I4/I3 local "
                "geometry but no rational section of height 79/12 is "
                "certified here."
            ),
            (
                "It is not yet identified with the X(6,79) rank-one seed "
                "or transported through the frozen neighbor chain."
            ),
            "This certificate does not improve the global rank-29 lower bound.",
        ],
        "implementation": {
            "language": "Python standard library",
            "script": "research/rank17_exact_split_seed_model.py",
        },
    }
    payload["certificate_sha256"] = hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    certificate = compute_certificate()
    if arguments.output:
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text())
        if committed != certificate:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
    print(
        json.dumps(
            {
                "status": certificate["truth_status"],
                "sha256": certificate["certificate_sha256"],
                "parameters": certificate["parameters"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
