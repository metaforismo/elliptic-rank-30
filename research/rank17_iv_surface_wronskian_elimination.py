#!/usr/bin/env sage -python
"""Eliminate the additive-IV surface curve using the Wronskian recurrence.

The sparse Wronskian bracket

    C=2AB+3(t-1)A'B-2(t-1)AB'

must have support only in degrees 3 and 4.  The equations at degrees 0,1,2
solve b1,b2,b3 from the low end.  Degrees 15,...,10 solve b9,...,b4 from the
high end.  The five remaining equations are coefficients 5,...,9, in the six
variables p0,a1,...,a5.  This script clears only denominators supported on
p0, saturates p0, and eliminates to the (a4,a5) plane.  Since

    s=a5-2,  r=a4-2*a5+1,

the resulting plane equation is also reported in (r,s).

This is exact modular algebra.  Boundary factors from C3*C4*B(1)=0 must be
classified after factorization.  No characteristic-zero or rank-30 claim is
made by one prime.
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


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def primitive_polynomial(value):
    polynomial = value
    if not polynomial:
        return polynomial
    leading = polynomial.lc()
    if leading != 1:
        polynomial = polynomial / leading
    return polynomial


def construct_reduced_equations(prime: int):
    field = GF(prime)
    base = PolynomialRing(
        field,
        names=("p0", "a1", "a2", "a3", "a4", "a5"),
        order="degrevlex",
    )
    p0, a1, a2, a3, a4, a5 = base.gens()
    fraction = base.fraction_field()
    time_ring = PolynomialRing(fraction, "t")
    t = time_ring.gen()
    A = time_ring(
        p0**2 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5 + t**6
    )

    coefficients = [fraction(p0**3)] + [None] * 9 + [fraction.one()]

    def build_B(target_index=None, target_value=None):
        values = []
        for index, coefficient in enumerate(coefficients):
            if index == target_index:
                values.append(fraction(target_value))
            elif coefficient is None:
                values.append(fraction.zero())
            else:
                values.append(coefficient)
        return time_ring(sum(values[index]*t**index for index in range(11)))

    def bracket(B):
        return 2*A*B + 3*(t-1)*A.derivative()*B - 2*(t-1)*A*B.derivative()

    def solve_linear(index: int, degree: int):
        zero_value = bracket(build_B(index, 0))[degree]
        one_value = bracket(build_B(index, 1))[degree]
        coefficient = one_value - zero_value
        if not coefficient:
            raise AssertionError(("vanishing recurrence coefficient", index, degree))
        solution = -zero_value / coefficient
        coefficients[index] = solution
        if bracket(build_B())[degree] != 0:
            raise AssertionError(("recurrence did not solve coefficient", index, degree))

    for index, degree in ((1, 0), (2, 1), (3, 2)):
        solve_linear(index, degree)
    for index, degree in ((9, 15), (8, 14), (7, 13), (6, 12), (5, 11), (4, 10)):
        solve_linear(index, degree)

    if any(coefficient is None for coefficient in coefficients):
        raise AssertionError("the bidirectional recurrence left an unknown B coefficient")
    B = build_B()
    C = bracket(B)
    forced_zero_degrees = [0, 1, 2, 10, 11, 12, 13, 14, 15]
    if any(C[degree] != 0 for degree in forced_zero_degrees):
        raise AssertionError("the recurrence replay lost a forced zero coefficient")

    residuals = []
    residual_metadata = []
    for degree in range(5, 10):
        numerator = C[degree].numerator()
        denominator = C[degree].denominator()
        polynomial = base(numerator)
        if not polynomial:
            raise AssertionError(("unexpected zero residual", degree))
        scalar = polynomial.lc()
        polynomial = base(polynomial / scalar)
        residuals.append(polynomial)
        residual_metadata.append({
            "wronskian_degree": degree,
            "total_degree": int(polynomial.total_degree()),
            "term_count": len(polynomial.dict()),
            "cleared_denominator": str(denominator),
            "normalized_leading_scalar": int(scalar),
            "sha256": hashlib.sha256(str(polynomial).encode()).hexdigest(),
        })

    open_data = {}
    for name, value in (
        ("C3", C[3]),
        ("C4", C[4]),
        ("B_at_one", B(1)),
    ):
        open_data[name] = {
            "numerator": str(value.numerator()),
            "denominator": str(value.denominator()),
        }

    coefficient_records = {
        f"b{index}": {
            "numerator": str(coefficients[index].numerator()),
            "denominator": str(coefficients[index].denominator()),
        }
        for index in range(1, 10)
    }
    return base, residuals, residual_metadata, coefficient_records, open_data


def factor_record(polynomial, first, second):
    return [
        {
            "factor": str(factor),
            "exponent": int(exponent),
            "total_degree": int(factor.total_degree()),
            "degree_first": int(factor.degree(first)),
            "degree_second": int(factor.degree(second)),
        }
        for factor, exponent in polynomial.factor()
    ]


def compute(prime: int, output: Path):
    base, residuals, residual_metadata, coefficient_records, open_data = (
        construct_reduced_equations(prime)
    )
    lex = PolynomialRing(
        GF(prime),
        names=("h", "p0", "a1", "a2", "a3", "a4", "a5"),
        order="lex",
    )
    h, p0, a1, a2, a3, a4, a5 = lex.gens()
    embedding = base.hom([p0, a1, a2, a3, a4, a5], lex)
    equations = [embedding(polynomial) for polynomial in residuals]
    equations.append(h*p0 - 1)

    initial = {
        "schema_version": 1,
        "status": "groebner_started",
        "truth_status": (
            f"exact Wronskian surface-curve elimination over F_{prime}; "
            "only p0 is saturated in the ideal; no characteristic-zero, "
            "section, or rank-30 conclusion"
        ),
        "prime": prime,
        "variable_order": [str(value) for value in lex.gens()],
        "equation_count": len(equations),
        "equation_degrees": [int(value.total_degree()) for value in equations],
        "equation_term_counts": [len(value.dict()) for value in equations],
        "residual_metadata": residual_metadata,
        "B_coefficient_recurrence": coefficient_records,
        "open_conditions_not_saturated": open_data,
        "sage_version": str(sage_version),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")

    ideal = lex.ideal(equations)
    std = singular_function("slimgb")
    dim_function = singular_function("dim")
    started = time.time()
    basis_raw = std(ideal)
    elapsed = time.time() - started
    dimension = int(dim_function(basis_raw))
    basis = [lex(value) for value in basis_raw]
    basis_text = "\n".join(str(value) for value in basis) + "\n"
    (output.parent / "wronskian-groebner-basis.txt").write_text(basis_text)

    keep = {a4, a5}
    pure = [
        polynomial for polynomial in basis
        if set(polynomial.variables()).issubset(keep)
    ]
    plane_ring = PolynomialRing(GF(prime), names=("a4", "a5"), order="lex")
    plane_polynomials = [plane_ring(str(value)) for value in pure]
    plane_basis = (
        list(plane_ring.ideal(plane_polynomials).groebner_basis())
        if plane_polynomials else []
    )
    plane_text = "\n".join(str(value) for value in plane_basis) + "\n"
    (output.parent / "wronskian-plane-basis-a4-a5.txt").write_text(plane_text)

    rs_ring = PolynomialRing(GF(prime), names=("r", "s"), order="lex")
    rr, ss = rs_ring.gens()
    transformed = [
        rs_ring(polynomial(rr + 2*ss + 3, ss + 2))
        for polynomial in plane_basis
    ]
    transformed_basis = (
        list(rs_ring.ideal(transformed).groebner_basis()) if transformed else []
    )
    rs_text = "\n".join(str(value) for value in transformed_basis) + "\n"
    (output.parent / "wronskian-plane-basis-r-s.txt").write_text(rs_text)

    result = dict(initial)
    result.update({
        "status": "completed",
        "elapsed_seconds": elapsed,
        "surface_krull_dimension": dimension,
        "groebner_basis_count": len(basis),
        "groebner_basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        "plane_basis_a4_a5": [
            {
                "polynomial": str(value),
                "total_degree": int(value.total_degree()),
                "degree_a4": int(value.degree(plane_ring.gen(0))),
                "degree_a5": int(value.degree(plane_ring.gen(1))),
                "factorization": factor_record(
                    value, plane_ring.gen(0), plane_ring.gen(1)
                ),
            }
            for value in plane_basis
        ],
        "plane_basis_r_s": [
            {
                "polynomial": str(value),
                "total_degree": int(value.total_degree()),
                "degree_r": int(value.degree(rr)),
                "degree_s": int(value.degree(ss)),
                "factorization": factor_record(value, rr, ss),
            }
            for value in transformed_basis
        ],
        "plane_basis_a4_a5_sha256": hashlib.sha256(plane_text.encode()).hexdigest(),
        "plane_basis_r_s_sha256": hashlib.sha256(rs_text.encode()).hexdigest(),
        "limitations": [
            "Only p0 is saturated during Groebner computation.",
            "Factors contained in C3*C4*B(1)=0 are boundary components and must be removed by exact component testing.",
            "A single modular plane equation is not an equation over Q.",
            "No height-79/12 section is imposed."
        ],
    })
    result["record_sha256"] = canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        result = compute(arguments.prime, arguments.output)
        print(json.dumps({
            key: result.get(key)
            for key in (
                "status", "prime", "surface_krull_dimension",
                "groebner_basis_count", "elapsed_seconds", "record_sha256"
            )
        }, indent=2, sort_keys=True), flush=True)
    except Exception as exc:
        error = {
            "schema_version": 1,
            "status": "error",
            "truth_status": "implementation or computation failure; no mathematical conclusion",
            "prime": arguments.prime,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
        }
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(error, indent=2, sort_keys=True) + "\n")
        print(json.dumps(error, indent=2, sort_keys=True), flush=True)
        raise


if __name__ == "__main__":
    main()
