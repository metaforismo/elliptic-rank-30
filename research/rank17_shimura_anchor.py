#!/usr/bin/env python3
"""Exact certificate for the Shimura-curve anchor behind the rank-17 K3.

Elkies identifies the moduli curve

    X(6,79)/<w_474>:  u^2 = 16 t^6 - 19 t^4 + 88 t^2 - 48

and singles out the rational orbit |t|=14/13, |u|=16064/2197 as
non-CM and as the source of an elliptic K3 surface of generic rank 17.

This script proves only the algebraic anchor that can be checked from the
published equation without reconstructing the moduli map:

* the sextic is squarefree, so the projective curve is nonsingular genus 2;
* the eight displayed affine rational points and two rational infinity points;
* the involution (t,u)->(-t,-u) has a nonsingular genus-1 quotient;
* exact images of the rational points on quartic and reciprocal cubic models.

It deliberately does NOT claim to prove that the |t|=14/13 orbit is non-CM,
that the associated K3 has Mordell-Weil rank 17, or that it specializes to
the rank-29 record. Those statements require the missing moduli/fibration map.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable

Poly = tuple[int, ...]  # coefficients in ascending order


def trim(poly: Iterable[int]) -> Poly:
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [0])


def derivative(poly: Poly) -> Poly:
    if len(poly) <= 1:
        return (0,)
    return trim(index * poly[index] for index in range(1, len(poly)))


def evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous_pivot = 1

    for column in range(size - 1):
        if work[column][column] == 0:
            swap_row = next(
                (
                    row
                    for row in range(column + 1, size)
                    if work[row][column] != 0
                ),
                None,
            )
            if swap_row is None:
                return 0
            work[column], work[swap_row] = work[swap_row], work[column]
            sign = -sign

        pivot = work[column][column]
        for row in range(column + 1, size):
            for target_column in range(column + 1, size):
                numerator = (
                    work[row][target_column] * pivot
                    - work[row][column] * work[column][target_column]
                )
                if numerator % previous_pivot:
                    raise AssertionError("Bareiss division was not exact")
                work[row][target_column] = numerator // previous_pivot
        previous_pivot = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
        for target_column in range(column + 1, size):
            work[column][target_column] = 0

    return sign * work[-1][-1]


def resultant(left: Poly, right: Poly) -> int:
    left = trim(left)
    right = trim(right)
    if left == (0,) or right == (0,):
        return 0

    degree_left = len(left) - 1
    degree_right = len(right) - 1
    size = degree_left + degree_right
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))

    sylvester: list[list[int]] = []
    for shift in range(degree_right):
        sylvester.append(
            [0] * shift
            + left_descending
            + [0] * (size - shift - len(left_descending))
        )
    for shift in range(degree_left):
        sylvester.append(
            [0] * shift
            + right_descending
            + [0] * (size - shift - len(right_descending))
        )
    return bareiss_determinant(sylvester)


def discriminant(poly: Poly) -> int:
    normalized = trim(poly)
    degree = len(normalized) - 1
    numerator = (
        (-1) ** (degree * (degree - 1) // 2)
        * resultant(normalized, derivative(normalized))
    )
    leading_coefficient = normalized[-1]
    if numerator % leading_coefficient:
        raise AssertionError("polynomial discriminant was not integral")
    return numerator // leading_coefficient


def fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


SEXTIC: Poly = (-48, 0, 88, 0, -19, 0, 16)
QUOTIENT_QUARTIC: Poly = (0, -48, 88, -19, 16)
RECIPROCAL_CUBIC: Poly = (16, -19, 88, -48)


def compute_certificate() -> dict[str, object]:
    sextic_resultant = resultant(SEXTIC, derivative(SEXTIC))
    sextic_discriminant = discriminant(SEXTIC)
    quartic_discriminant = discriminant(QUOTIENT_QUARTIC)
    cubic_discriminant = discriminant(RECIPROCAL_CUBIC)

    if sextic_resultant != -960467703986795839488:
        raise AssertionError("unexpected sextic resultant")
    if sextic_discriminant != 60029231499174739968:
        raise AssertionError("unexpected sextic discriminant")
    if quartic_discriminant != -80518053888:
        raise AssertionError("unexpected quotient-quartic discriminant")
    if cubic_discriminant != -34947072:
        raise AssertionError("unexpected reciprocal-cubic discriminant")

    positive_representatives = [
        (Fraction(2), Fraction(32)),
        (Fraction(14, 13), Fraction(16064, 2197)),
    ]
    affine_points: list[list[str]] = []
    quartic_images: dict[tuple[str, str], list[str]] = {}
    cubic_images: dict[tuple[str, str], list[str]] = {}

    for t_absolute, u_absolute in positive_representatives:
        if u_absolute * u_absolute != evaluate(SEXTIC, t_absolute):
            raise AssertionError("distinguished rational point is off the curve")

        for t_sign in (1, -1):
            for u_sign in (1, -1):
                t_value = t_sign * t_absolute
                u_value = u_sign * u_absolute
                if u_value * u_value != evaluate(SEXTIC, t_value):
                    raise AssertionError("rational orbit point is off the curve")

                source_key = (
                    fraction_string(t_value),
                    fraction_string(u_value),
                )
                affine_points.append(list(source_key))

                quotient_x = t_value * t_value
                quotient_y = t_value * u_value
                if quotient_y * quotient_y != evaluate(
                    QUOTIENT_QUARTIC, quotient_x
                ):
                    raise AssertionError("quartic quotient image is invalid")
                quartic_images[source_key] = [
                    fraction_string(quotient_x),
                    fraction_string(quotient_y),
                ]

                reciprocal_x = 1 / quotient_x
                reciprocal_y = u_value / (t_value**3)
                if reciprocal_y * reciprocal_y != evaluate(
                    RECIPROCAL_CUBIC, reciprocal_x
                ):
                    raise AssertionError("reciprocal cubic image is invalid")
                cubic_images[source_key] = [
                    fraction_string(reciprocal_x),
                    fraction_string(reciprocal_y),
                ]

    if len({tuple(point) for point in affine_points}) != 8:
        raise AssertionError("expected eight distinct affine rational points")

    # The sextic is even, so iota(t,u)=(-t,-u) preserves the equation.
    for t_text, _u_text in affine_points:
        t_value = Fraction(t_text)
        if evaluate(SEXTIC, -t_value) != evaluate(SEXTIC, t_value):
            raise AssertionError("the proposed involution does not preserve C")

    exact_checks = {
        "F(2)=32^2": evaluate(SEXTIC, Fraction(2)) == Fraction(32**2),
        "F(14/13)=(16064/2197)^2": (
            evaluate(SEXTIC, Fraction(14, 13))
            == Fraction(16064, 2197) ** 2
        ),
        "N=6*79": 6 * 79 == 474,
    }
    if not all(exact_checks.values()):
        raise AssertionError(f"distinguished checks failed: {exact_checks}")

    sorted_affine_points = sorted(
        affine_points,
        key=lambda point: (Fraction(point[0]), Fraction(point[1])),
    )

    payload: dict[str, object] = {
        "schema_version": 1,
        "anchor_id": "elkies_rank17_shimura_anchor_N474",
        "exact_claim": (
            "The genus-2 curve C: u^2=16t^6-19t^4+88t^2-48 is "
            "nonsingular, has the listed rational affine points and two "
            "rational points at infinity, and the involution "
            "(t,u)->(-t,-u) has a nonsingular genus-1 quotient."
        ),
        "shimura_metadata": {
            "N": 474,
            "factorization": [6, 79],
            "curve_label": "X(6,79)/<w_474>",
            "status": (
                "identification from Elkies's primary-source lectures/paper; "
                "the moduli map is not reconstructed by this certificate"
            ),
        },
        "curve": {
            "equation": "u^2=16*t^6-19*t^4+88*t^2-48",
            "polynomial_coefficients_ascending": [
                str(coefficient) for coefficient in SEXTIC
            ],
            "degree": 6,
            "polynomial_resultant_with_derivative": str(sextic_resultant),
            "polynomial_discriminant": str(sextic_discriminant),
            "squarefree": sextic_discriminant != 0,
            "genus": 2,
            "rational_points_at_infinity_weighted_projective_T_Z_U": [
                ["1", "0", "4"],
                ["1", "0", "-4"],
            ],
        },
        "rational_affine_points": sorted_affine_points,
        "distinguished_exact_checks": exact_checks,
        "bielliptic_involution": {
            "map": "iota(t,u)=(-t,-u)",
            "order": 2,
            "invariants": {"x": "t^2", "y": "t*u"},
            "quotient_equation": (
                "y^2=16*x^4-19*x^3+88*x^2-48*x"
            ),
            "quotient_polynomial_coefficients_ascending": [
                str(coefficient) for coefficient in QUOTIENT_QUARTIC
            ],
            "quotient_discriminant": str(quartic_discriminant),
            "quotient_nonsingular": quartic_discriminant != 0,
            "quotient_genus": 1,
            "reciprocal_cubic_model": (
                "Y^2=-48*X^3+88*X^2-19*X+16, "
                "X=1/t^2, Y=u/t^3"
            ),
            "reciprocal_cubic_discriminant": str(cubic_discriminant),
        },
        "point_images": {
            "quartic_quotient_by_source_point": {
                ",".join(key): value
                for key, value in sorted(quartic_images.items())
            },
            "reciprocal_cubic_by_source_point": {
                ",".join(key): value
                for key, value in sorted(cubic_images.items())
            },
        },
        "source_claims_not_proved_here": [
            "The |t|=14/13 orbit is non-CM.",
            (
                "That non-CM moduli point yields an elliptic K3 surface "
                "over Q(s) with Mordell-Weil rank 17."
            ),
            (
                "This rank-17 surface is the geometric source used in the "
                "rank-28 and rank-29 record searches."
            ),
        ],
        "missing_bridge": (
            "An explicit moduli map from this Shimura point to a Weierstrass "
            "K3 fibration, the 17 generic sections, and the specialization "
            "map to the published rank-29 curve."
        ),
        "conditional_assumptions": [],
        "implementation": {
            "language": "Python standard library",
            "script": "research/rank17_shimura_anchor.py",
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the freshly recomputed JSON certificate",
    )
    parser.add_argument(
        "--compare",
        type=Path,
        help="compare recomputation with a committed JSON certificate",
    )
    arguments = parser.parse_args()

    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != certificate:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
        print(f"matched {arguments.compare}")

    print("N=474=6*79")
    print("Shimura anchor polynomial: squarefree degree 6, genus 2")
    print("rational affine points verified: 8; rational infinity points: 2")
    print("bielliptic quotient: nonsingular genus 1")
    print(f"certificate sha256: {certificate['certificate_sha256']}")
    print("NOT proved here: non-CM status or rank-17 K3 moduli map")


if __name__ == "__main__":
    main()
