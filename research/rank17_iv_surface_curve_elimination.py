#!/usr/bin/env python3
"""Eliminate the normalized I12+I4+IV surface locus to a plane curve.

The auxiliary jet-sign involution lets us take e0=1 without changing c4 or c6.
For a requested finite field and one of the invariant projections

    (r,s), (p0,r), (p0,s),

this script constructs the five exact surface equations, saturates p0 and the
coefficients enforcing exact I4, IV, and I12 orders, computes a lexicographic
Groebner basis, and extracts the corresponding plane elimination ideal.

This is exact modular algebra.  A stable plane equation across several good
primes can later be reconstructed over Q, but no such reconstruction is claimed
by a single run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import traceback
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version

PROJECTIONS = {
    "r-s": ("r", "s"),
    "p0-r": ("p0", "r"),
    "p0-s": ("p0", "s"),
}


def convolution(left, right, length=None):
    if length is None:
        length = len(left) + len(right) - 1
    zero = left[0] * 0
    result = [zero for _ in range(length)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            target = left_index + right_index
            if target >= length:
                break
            result[target] += left_value * right_value
    return result


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def factor_record(polynomial) -> list[dict[str, object]]:
    return [
        {
            "factor": str(factor),
            "total_degree": int(factor.total_degree()),
            "degree_first": int(factor.degree(factor.parent().gen(0))),
            "degree_second": int(factor.degree(factor.parent().gen(1))),
            "exponent": int(exponent),
        }
        for factor, exponent in polynomial.factor()
    ]


def construct_surface_equations(ring):
    variables = {str(generator): generator for generator in ring.gens()}
    p0 = variables["p0"]
    p1 = variables["p1"]
    p2 = variables["p2"]
    p3 = variables["p3"]
    r = variables["r"]
    s = variables["s"]
    h = variables["h"]

    a0 = p0**2
    a1 = 2*p0*p1
    a2 = 2*p0*p2 + p1**2
    a3 = 2*p0*p3 + 2*p1*p2
    l0 = a0 + a1 + a2 + a3 + s + 1
    l1 = a1 + 2*a2 + 3*a3 + 7*s + 8
    a4 = l1 - 5*l0 + r
    a5 = 4*l0 - l1 - 2*r
    c4 = [a0, a1, a2, a3, a4, a5, r, s, ring.one()]

    reversed_c4 = list(reversed(c4))
    square = convolution(reversed_c4, reversed_c4, 12)
    cube = convolution(square, reversed_c4, 12)
    inverse_two = ring.base_ring()(2) ** -1
    root = [ring.zero() for _ in range(12)]
    root[0] = ring.one()
    for order in range(1, 12):
        correction = sum(
            (root[index] * root[order-index] for index in range(1, order)),
            ring.zero(),
        )
        root[order] = (cube[order] - correction) * inverse_two

    c6 = [ring.zero() for _ in range(13)]
    c6[0] = p0**3
    c6[12] = ring.one()
    for order in range(1, 12):
        c6[12-order] = root[order]

    equations = [
        c6[1] - 3*p0**2*p1,
        c6[2] - 3*(p0**2*p2 + p0*p1**2),
        c6[3] - (3*p0**2*p3 + 6*p0*p1*p2 + p1**3),
        sum(c6, ring.zero()),
        sum((index*c6[index] for index in range(13)), ring.zero()),
    ]

    c4_square = convolution(c4, c4)
    c4_cube = convolution(c4_square, c4)
    c6_square = convolution(c6, c6)
    delta = [
        c4_cube[index] - c6_square[index]
        for index in range(25)
    ]
    b_at_one = sum(
        (
            (index*(index-1)//2) * c6[index]
            for index in range(2, 13)
        ),
        ring.zero(),
    )
    open_product = p0 * b_at_one * delta[4] * delta[12]
    equations.append(h * open_product - 1)
    return {
        "equations": equations,
        "c4": c4,
        "c6": c6,
        "delta": delta,
        "b_at_one": b_at_one,
        "open_product": open_product,
    }


def compute(prime: int, projection_name: str, output: Path) -> dict[str, object]:
    if projection_name not in PROJECTIONS:
        raise ValueError(projection_name)
    keep = PROJECTIONS[projection_name]
    eliminate = [
        name for name in ("p0", "p1", "p2", "p3", "r", "s", "h")
        if name not in keep
    ]
    variable_names = tuple(eliminate + list(keep))
    field = GF(prime)
    ring = PolynomialRing(field, names=variable_names, order="lex")
    construction = construct_surface_equations(ring)
    equations = construction["equations"]

    initial = {
        "status": "groebner_started",
        "truth_status": (
            f"exact saturated surface-curve elimination over F_{prime}; "
            "no characteristic-zero, section, or rank-30 conclusion"
        ),
        "prime": prime,
        "projection": projection_name,
        "kept_coordinates": list(keep),
        "variable_order": list(variable_names),
        "equation_count": len(equations),
        "equation_degrees": [int(equation.total_degree()) for equation in equations],
        "sage_version": str(sage_version),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")

    started = time.time()
    ideal = ring.ideal(equations)
    std = singular_function("std")
    dim_function = singular_function("dim")
    basis_raw = std(ideal)
    elapsed = time.time() - started
    basis = [ring(polynomial) for polynomial in basis_raw]
    dimension = int(dim_function(basis_raw))
    basis_text = "\n".join(str(polynomial) for polynomial in basis) + "\n"
    (output.parent / "groebner-basis.txt").write_text(basis_text)

    keep_generators = {ring(name) for name in keep}
    pure = [
        polynomial for polynomial in basis
        if set(polynomial.variables()).issubset(keep_generators)
    ]
    pure_text = "\n".join(str(polynomial) for polynomial in pure) + "\n"
    (output.parent / "plane-elimination-polynomials.txt").write_text(pure_text)

    plane_ring = PolynomialRing(field, names=keep, order="lex")
    plane_polynomials = [
        plane_ring(str(polynomial)) for polynomial in pure
    ]
    plane_basis = (
        list(plane_ring.ideal(plane_polynomials).groebner_basis())
        if plane_polynomials else []
    )
    plane_basis_text = "\n".join(str(polynomial) for polynomial in plane_basis) + "\n"
    (output.parent / "plane-groebner-basis.txt").write_text(plane_basis_text)

    result = dict(initial)
    result.update({
        "status": "completed",
        "elapsed_seconds": elapsed,
        "surface_krull_dimension": dimension,
        "groebner_basis_count": len(basis),
        "groebner_basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        "pure_plane_polynomial_count": len(pure),
        "plane_groebner_basis_count": len(plane_basis),
        "plane_groebner_basis_sha256": hashlib.sha256(
            plane_basis_text.encode()
        ).hexdigest(),
        "plane_basis": [
            {
                "polynomial": str(polynomial),
                "total_degree": int(polynomial.total_degree()),
                "degree_first": int(polynomial.degree(plane_ring.gen(0))),
                "degree_second": int(polynomial.degree(plane_ring.gen(1))),
                "factorization": factor_record(polynomial),
            }
            for polynomial in plane_basis
        ],
    })
    result["record_sha256"] = canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--projection", choices=sorted(PROJECTIONS), required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    try:
        result = compute(arguments.prime, arguments.projection, arguments.output)
        print(json.dumps({
            key: result.get(key)
            for key in (
                "status",
                "prime",
                "projection",
                "surface_krull_dimension",
                "plane_groebner_basis_count",
                "elapsed_seconds",
                "record_sha256",
            )
        }, indent=2, sort_keys=True), flush=True)
    except Exception as exc:
        error = {
            "status": "error",
            "truth_status": "implementation or computation failure; no mathematical conclusion",
            "error": repr(exc),
            "traceback": traceback.format_exc(),
        }
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(error, indent=2, sort_keys=True) + "\n")
        print(json.dumps(error, indent=2, sort_keys=True), flush=True)
        raise


if __name__ == "__main__":
    main()
