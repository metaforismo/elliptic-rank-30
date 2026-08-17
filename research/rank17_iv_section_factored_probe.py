#!/usr/bin/env sage -python
"""Factored exact section-incidence probe on a split IV surface over F_11.

Use

    c4=(t-1)^2*A,  c6=(t-1)^2*B,
    X=(t-1)*W,     Y=t*(t-1)*Z.

The section identity then reduces exactly to

    Z^2 = ((t-1)*(W^3-3*A*W*D^4)-2*B*D^6)/t^2.

The numerator is divisible by t^2 once W(0) is fixed at the I4 node.  A
square-root recurrence from t=0 and t=infinity eliminates z1,...,z9 and leaves
an overdetermined ideal in ten variables.  This is a finite-field probe only.
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


def divide_by_t_minus_one_squared(coefficients, ring):
    remainder = list(coefficients)
    quotient = [ring(0) for _ in range(len(coefficients) - 2)]
    # (t-1)^2 has ascending coefficients [1,-2,1].
    for degree in range(len(coefficients) - 1, 1, -1):
        coefficient = remainder[degree]
        quotient[degree - 2] = coefficient
        remainder[degree - 2] -= coefficient
        remainder[degree - 1] += 2 * coefficient
        remainder[degree] -= coefficient
    if remainder[0] != 0 or remainder[1] != 0:
        raise AssertionError("surface invariant is not divisible by (t-1)^2")
    return quotient


def build_system(surface_index: int, q_value: int, i4_sign: int, iv_sign: int):
    prime = 11
    field = GF(prime)
    surface = SURFACES_P11[surface_index]
    names = ("d0", "d1", "w1", "w2", "w3", "w4", "w5", "w6", "h", "k")
    ring = PolynomialRing(field, names=names, order="degrevlex")
    d0, d1, w1, w2, w3, w4, w5, w6, h, k = ring.gens()

    q = field(q_value)
    if q == 0:
        raise ValueError("the first factored probe uses the y12 != 0 infinity chart")
    w7 = ring(q**2 + 2)
    z10 = ring(q * (q**2 + 3))
    if z10 == 0:
        raise ValueError("the selected infinity point has zero y-coordinate")

    p0, p1 = map(field, surface["parameters"][:2])
    e0 = field(surface["e0"])
    alpha = field(i4_sign * surface["i4_root"])
    beta = field(iv_sign * surface["iv_root"])
    c4 = [ring(field(value)) for value in surface["c4"]]
    c6 = [ring(field(value)) for value in surface["c6"]]
    A = divide_by_t_minus_one_squared(c4, ring)
    B = divide_by_t_minus_one_squared(c6, ring)

    denominator = [d0, d1, ring(1)]
    denominator_squared = convolution(denominator, denominator)
    denominator_fourth = convolution(denominator_squared, denominator_squared)
    denominator_sixth = convolution(denominator_fourth, denominator_squared)

    # X=(t-1)W and X(0)=-e0*p0*D(0)^2 force W(0)=e0*p0*d0^2.
    w0 = e0 * p0 * d0**2
    W = [w0, w1, w2, w3, w4, w5, w6, w7]
    w_squared = convolution(W, W)
    w_cubed = convolution(w_squared, W)
    A_W = convolution(A, W)
    A_W_D4 = convolution(A_W, denominator_fourth)
    core = [w_cubed[index] - 3 * A_W_D4[index] for index in range(22)]
    times_t_minus_one = [ring(0) for _ in range(23)]
    for index, coefficient in enumerate(core):
        times_t_minus_one[index] -= coefficient
        times_t_minus_one[index + 1] += coefficient
    B_D6 = convolution(B, denominator_sixth)
    bracket = [
        times_t_minus_one[index] - 2 * B_D6[index]
        for index in range(23)
    ]
    if bracket[0] != 0 or bracket[1] != 0:
        raise AssertionError("I4 node substitution did not produce the t^2 factor")
    square_target = bracket[2:]
    if len(square_target) != 21:
        raise AssertionError("unexpected factored square-target degree")

    # Y=t(t-1)Z, so Y_1=-Z_0.  The I4 tangent condition therefore fixes z0.
    x1 = w0 - w1
    moving_node_first_jet = x1 + e0 * (
        p0 * denominator_squared[1] + p1 * denominator_squared[0]
    )
    z = [None for _ in range(11)]
    z[0] = -alpha * d0 * moving_node_first_jet
    z[10] = z10

    if square_target[0] - z[0] ** 2 != 0:
        raise AssertionError("the I4 tangent square did not cancel identically")
    if square_target[20] - z[10] ** 2 != 0:
        raise AssertionError("the infinity tangent square did not cancel identically")

    # h=1/(2*z0) simultaneously removes d0=0 and the identity tangent z0=0.
    for degree in range(1, 6):
        known = ring(0)
        for index in range(1, degree):
            known += z[index] * z[degree - index]
        z[degree] = h * (square_target[degree] - known)

    inverse_two_z10 = (2 * field(z10)) ** -1
    overlap_equation = None
    for degree in range(19, 14, -1):
        target_index = degree - 10
        known = ring(0)
        for index in range(target_index + 1, 10):
            known += z[index] * z[degree - index]
        high_value = inverse_two_z10 * (square_target[degree] - known)
        if z[target_index] is None:
            z[target_index] = high_value
        else:
            overlap_equation = z[target_index] - high_value
    if overlap_equation is None or any(value is None for value in z):
        raise AssertionError("factored two-sided recurrence is incomplete")

    d_at_one = d0 + d1 + 1
    equations = [2 * z[0] * h - 1, d_at_one * k - 1, overlap_equation]
    for degree in range(6, 15):
        square_coefficient = ring(0)
        for index in range(11):
            other = degree - index
            if 0 <= other < 11:
                square_coefficient += z[index] * z[other]
        equations.append(square_coefficient - square_target[degree])

    # At t=1, Y'(1)=Z(1); this selects one of the two outer IV components.
    equations.append(sum(z, ring(0)) - beta * d_at_one**3)

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
        "automatic_zero_target_degrees": [0, 20],
        "middle_checked_degrees": list(range(6, 15)),
        "factored_identity": (
            "Z^2=((t-1)*(W^3-3*A*W*D^4)-2*B*D^6)/t^2"
        ),
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
            "EXACT factored finite-field section-incidence probe; "
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
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    if result["status"] == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
