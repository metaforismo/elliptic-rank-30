#!/usr/bin/env sage -python
"""Reduced exact section-incidence probe on a split IV surface over F_p.

For a fixed normalized surface, denominator chart, smooth point at infinity,
and local tangent signs, the polynomial-square identity is solved from both
ends.  The high and low square-root recurrences eliminate y_2,...,y_11,
leaving a small overdetermined ideal in the denominator and x-coefficients.

This is a finite-field research probe.  It makes no characteristic-zero or
rank-30 claim.
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

from sage.all import GF, PolynomialRing, version

SURFACES_P11 = [
    {
        "parameters": [1, 3, 1, 4, 9, 5],
        "e0": -1,
        "c4": [1, 6, 0, 3, 9, 10, 9, 5, 1],
        "c6": [10, 2, 3, 9, 5, 2, 2, 8, 7, 10, 5, 2, 1],
        "i4_root": 5,
        "iv_root": 4,
    },
    {
        "parameters": [2, 3, 7, 4, 0, 1],
        "e0": 1,
        "c4": [4, 1, 4, 3, 7, 1, 0, 1, 1],
        "c6": [8, 3, 6, 8, 1, 6, 1, 1, 5, 9, 10, 7, 1],
        "i4_root": 4,
        "iv_root": 3,
    },
    {
        "parameters": [5, 6, 0, 5, 7, 2],
        "e0": -1,
        "c4": [3, 5, 3, 6, 2, 4, 7, 2, 1],
        "c6": [7, 1, 10, 3, 9, 4, 5, 0, 6, 5, 1, 3, 1],
        "i4_root": 2,
        "iv_root": 2,
    },
    {
        "parameters": [5, 9, 1, 8, 5, 6],
        "e0": -1,
        "c4": [3, 2, 3, 10, 1, 2, 5, 6, 1],
        "c6": [7, 7, 8, 7, 8, 4, 3, 1, 0, 1, 10, 9, 1],
        "i4_root": 2,
        "iv_root": 2,
    },
]


def convolution(left, right):
    output = [left[0] * 0 for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i + j] += a * b
    return output


def build_system(surface_index: int, q_value: int, i4_sign: int, iv_sign: int):
    prime = 11
    field = GF(prime)
    surface = SURFACES_P11[surface_index]
    names = ("d0", "d1", "x2", "x3", "x4", "x5", "x6", "x7", "h", "k")
    ring = PolynomialRing(field, names=names, order="degrevlex")
    d0, d1, x2, x3, x4, x5, x6, x7, h, k = ring.gens()

    q = field(q_value)
    if q == 0:
        raise ValueError("this first probe requires nonzero y-leading coefficient")
    x8 = ring(q**2 + 2)
    y12 = ring(q * (q**2 + 3))
    if y12 == 0:
        raise ValueError("the selected infinity chart has y12=0")

    p0, p1 = map(field, surface["parameters"][:2])
    e0 = field(surface["e0"])
    alpha = field(i4_sign * surface["i4_root"])
    beta = field(iv_sign * surface["iv_root"])
    c4 = [ring(field(value)) for value in surface["c4"]]
    c6 = [ring(field(value)) for value in surface["c6"]]

    denominator = [d0, d1, ring(1)]
    denominator_squared = convolution(denominator, denominator)
    denominator_fourth = convolution(denominator_squared, denominator_squared)
    denominator_sixth = convolution(denominator_fourth, denominator_squared)

    singular_x0 = -e0 * p0
    x0 = singular_x0 * d0**2
    x1 = -(x0 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
    x_polynomial = [x0, x1, x2, x3, x4, x5, x6, x7, x8]

    x_squared = convolution(x_polynomial, x_polynomial)
    x_cubed = convolution(x_squared, x_polynomial)
    c4_x = convolution(c4, x_polynomial)
    c4_x_d4 = convolution(c4_x, denominator_fourth)
    c6_d6 = convolution(c6, denominator_sixth)
    rhs = [
        x_cubed[index] - 3 * c4_x_d4[index] - 2 * c6_d6[index]
        for index in range(25)
    ]

    moving_node_first_jet = x1 + e0 * (
        p0 * denominator_squared[1] + p1 * denominator_squared[0]
    )
    y = [None for _ in range(13)]
    y[0] = ring(0)
    y[1] = alpha * d0 * moving_node_first_jet
    y[12] = y12

    automatic_zero_checks = [
        rhs[0],
        rhs[1],
        rhs[2] - y[1] ** 2,
        rhs[24] - y[12] ** 2,
    ]
    if any(value != 0 for value in automatic_zero_checks):
        raise AssertionError("local or infinity coefficient identities did not cancel")

    # The inverse h=1/(2*y1) makes the low recurrence polynomial and also
    # saturates away d0=0 and the identity-component tangent u0=0.
    for degree in range(3, 8):
        known = ring(0)
        for index in range(2, degree - 1):
            known += y[index] * y[degree - index]
        y[degree - 1] = h * (rhs[degree] - known)

    inverse_two_y12 = (2 * field(y12)) ** -1
    overlap_equation = None
    for degree in range(23, 17, -1):
        target_index = degree - 12
        known = ring(0)
        for index in range(target_index + 1, 12):
            known += y[index] * y[degree - index]
        high_value = inverse_two_y12 * (rhs[degree] - known)
        if y[target_index] is None:
            y[target_index] = high_value
        else:
            overlap_equation = y[target_index] - high_value
    if overlap_equation is None or any(value is None for value in y):
        raise AssertionError("two-sided square-root recurrence is incomplete")

    equations = [2 * y[1] * h - 1, (d0 + d1 + 1) * k - 1, overlap_equation]
    for degree in range(8, 18):
        square_coefficient = ring(0)
        for index in range(13):
            other = degree - index
            if 0 <= other < 13:
                square_coefficient += y[index] * y[other]
        equations.append(square_coefficient - rhs[degree])

    # These two equations select an actual outer IV component rather than only
    # retaining the nonreduced implication Y(1)^2=0 in the quotient algebra.
    y_at_one = sum(y, ring(0))
    y_derivative_at_one = sum((index * y[index] for index in range(13)), ring(0))
    d_at_one = d0 + d1 + 1
    equations.extend([
        y_at_one,
        y_derivative_at_one - beta * d_at_one**3,
    ])

    metadata = {
        "prime": prime,
        "surface_index": surface_index,
        "surface": surface,
        "q_at_infinity": q_value,
        "i4_tangent_sign": i4_sign,
        "iv_tangent_sign": iv_sign,
        "variable_names": list(names),
        "equation_count": len(equations),
        "equation_degrees": [int(polynomial.total_degree()) for polynomial in equations],
        "equation_term_counts": [len(polynomial.dict()) for polynomial in equations],
        "automatic_zero_coefficient_degrees": [0, 1, 2, 24],
        "middle_checked_degrees": list(range(8, 18)),
    }
    return ring, equations, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", type=int, default=0)
    parser.add_argument("--q", type=int, default=1)
    parser.add_argument("--i4-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--iv-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    result = {
        "schema_version": 1,
        "truth_status": (
            "EXACT finite-field reduced section-incidence probe; "
            "no characteristic-zero or rank-30 conclusion"
        ),
        "sage_version": str(version()),
        "status": "started",
    }
    try:
        ring, equations, metadata = build_system(
            arguments.surface,
            arguments.q,
            arguments.i4_sign,
            arguments.iv_sign,
        )
        result.update(metadata)
        ideal = ring.ideal(equations)
        started = time.time()
        basis = ideal.groebner_basis()
        elapsed = time.time() - started
        dimension = int(ideal.dimension())
        basis_text = "\n".join(str(polynomial) for polynomial in basis) + "\n"
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        (arguments.output.parent / "basis.txt").write_text(basis_text)
        result.update({
            "status": "completed",
            "elapsed_seconds": elapsed,
            "ideal_dimension": dimension,
            "basis_count": len(basis),
            "basis_degrees": [int(polynomial.total_degree()) for polynomial in basis],
            "basis_term_counts": [len(polynomial.dict()) for polynomial in basis],
            "unit_ideal": any(polynomial == 1 for polynomial in basis),
        })
        if dimension == 0 and not result["unit_ideal"]:
            try:
                result["quotient_vector_space_dimension"] = int(
                    ideal.vector_space_dimension()
                )
            except Exception as exc:
                result["quotient_dimension_error"] = repr(exc)
            try:
                points = ideal.variety(ring=GF(11))
                result["rational_point_count"] = len(points)
                result["rational_points"] = [
                    {str(key): int(value) for key, value in point.items()}
                    for point in points
                ]
            except Exception as exc:
                result["variety_error"] = repr(exc)
    except Exception as exc:
        result.update({
            "status": "error",
            "error": repr(exc),
            "traceback": traceback.format_exc(),
        })
    arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    if result["status"] == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
