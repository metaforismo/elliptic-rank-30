#!/usr/bin/env sage -python
"""Eliminate the open additive-IV surface curve in three reduced variables.

This wrapper reuses the exact logarithmic-derivative residual construction and
adds the missing open condition which separates genuine surfaces from the
boundary component with

    (t-1)^2 A^3 - B^2 = 0.

Before division by K, put

    J=(t-1)^2*N^3*K^2-T^2.

On the differential locus with A=N/K^2 and B=T/K^4 one has

    J = kappa*t^3*P*K^8.

Because P is monic and K=-12*t+k, the t^16 coefficient is
kappa*(-12)^8.  Its numerator therefore gives an exact open scalar.  The ideal
is saturated by both the lower-recurrence determinant and this scalar before
lexicographic elimination to the (r2,r3) plane.

A modular result is not an equation over Q.  No section or rank-30 conclusion
is made.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
import traceback
from pathlib import Path

from sage.all import PolynomialRing
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "rank17_iv_three_variable_elimination_sage.py"


def load_implementation():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_three_variable_elimination_sage", IMPLEMENTATION
    )
    if spec is None or spec.loader is None:
        raise ImportError(IMPLEMENTATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def construct_open_scalar(construction):
    base = construction["base"]
    r2, r3, k = base.gens()
    fraction = base.fraction_field()
    time_ring = PolynomialRing(fraction, "t")
    t = time_ring.gen()
    s = t - 1
    lower = construction["lower_solution"]
    m0 = fraction(lower["m0"])
    r1 = fraction(lower["r1"])
    r0 = fraction(lower["r0"])
    m1 = fraction(lower["m1"])
    m2 = fraction(lower["m2"])

    R = r0 + r1*t + r2*t**2 + r3*t**3 + t**4
    P = t*R
    Q = P.derivative() + 3*R
    K = k - 12*t
    M = m0 + m1*t + m2*t**2 + 80*t**3
    N = Q**2 + M*P
    T = (
        N*K*(s*Q - 2*P)
        - 3*P*s*(N.derivative()*K - 2*N*K.derivative())
    )
    J = s**2*N**3*K**2 - T**2
    coefficient = J[16] / ((-12)**8)
    numerator = base(coefficient.numerator())
    denominator = base(coefficient.denominator())
    if not numerator:
        raise AssertionError("the discriminant scalar vanished identically")

    determinant = construction["determinant"]
    removed_power = 0
    while True:
        quotient, remainder = numerator.quo_rem(determinant)
        if remainder:
            break
        numerator = quotient
        removed_power += 1
    numerator = numerator.monic()
    return {
        "polynomial": numerator,
        "denominator": denominator,
        "removed_determinant_power": removed_power,
        "total_degree": int(numerator.total_degree()),
        "degree_r2": int(numerator.degree(r2)),
        "degree_r3": int(numerator.degree(r3)),
        "degree_k": int(numerator.degree(k)),
        "term_count": len(numerator.dict()),
        "sha256": hashlib.sha256(str(numerator).encode()).hexdigest(),
        "identity": "coeff_t16((t-1)^2*N^3*K^2-T^2)=kappa*(-12)^8 on the differential locus",
    }


def replay_open_source(implementation, census_path, construction, open_scalar):
    replay = implementation.replay_source(census_path, construction)
    base = construction["base"]
    r2, r3, k = base.gens()
    values = []
    for point in replay["points"]:
        substitution = {
            r2: point["r2"],
            r3: point["r3"],
            k: point["k"],
        }
        value = int(open_scalar["polynomial"](substitution))
        if value == 0:
            raise AssertionError((
                "certified source surface lies on the zero-discriminant boundary",
                point,
            ))
        values.append(value)
    replay["all_source_points_avoid_zero_discriminant_boundary"] = True
    replay["open_scalar_values"] = values
    return replay


def compute(prime: int, census_path: Path, output: Path):
    implementation = load_implementation()
    construction = implementation.construct_residuals(prime)
    open_scalar = construct_open_scalar(construction)
    source_replay = replay_open_source(
        implementation, census_path, construction, open_scalar
    )
    field = construction["field"]
    base = construction["base"]

    lex = PolynomialRing(
        field, names=("h", "k", "r2", "r3"), order="lex"
    )
    h, k, r2, r3 = lex.gens()
    embedding = base.hom([r2, r3, k], lex)
    equations = [embedding(value) for value in construction["residuals"]]
    saturation_product = (
        embedding(construction["determinant"])
        * embedding(open_scalar["polynomial"])
    )
    equations.append(h*saturation_product - 1)

    initial = {
        "schema_version": 1,
        "status": "groebner_started",
        "truth_status": (
            f"exact determinant- and discriminant-saturated three-variable "
            f"additive-IV surface elimination over F_{prime}; no characteristic-zero, "
            "section, or rank-30 conclusion"
        ),
        "prime": prime,
        "sage_version": str(sage_version),
        "variable_order": [str(value) for value in lex.gens()],
        "residual_metadata": construction["residual_metadata"],
        "determinant": str(construction["determinant"]),
        "determinant_ratio": construction["determinant_ratio"],
        "open_discriminant_scalar": {
            key: value
            for key, value in open_scalar.items()
            if key != "polynomial" and key != "denominator"
        } | {
            "polynomial": str(open_scalar["polynomial"]),
            "denominator": str(open_scalar["denominator"]),
        },
        "source_replay": source_replay,
        "equation_count": len(equations),
        "equation_degrees": [int(value.total_degree()) for value in equations],
        "equation_term_counts": [len(value.dict()) for value in equations],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")

    ideal = lex.ideal(equations)
    slimgb = singular_function("slimgb")
    dim_function = singular_function("dim")
    started = time.time()
    basis_raw = slimgb(ideal)
    elapsed = time.time() - started
    dimension = int(dim_function(basis_raw))
    basis = [lex(value) for value in basis_raw]
    basis_text = "\n".join(str(value) for value in basis) + "\n"
    (output.parent / "three-variable-open-lex-basis.txt").write_text(basis_text)

    keep = {r2, r3}
    pure = [
        polynomial for polynomial in basis
        if set(polynomial.variables()).issubset(keep)
    ]
    plane = PolynomialRing(field, names=("r2", "r3"), order="lex")
    x, y = plane.gens()
    plane_polynomials = [plane(str(value)) for value in pure]
    plane_basis = (
        list(plane.ideal(plane_polynomials).groebner_basis())
        if plane_polynomials else []
    )
    plane_text = "\n".join(str(value) for value in plane_basis) + "\n"
    (output.parent / "three-variable-open-plane-r2-r3.txt").write_text(plane_text)

    result = dict(initial)
    result.update({
        "status": "completed",
        "elapsed_seconds": elapsed,
        "saturated_krull_dimension": dimension,
        "groebner_basis_count": len(basis),
        "groebner_basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        "plane_basis": [
            {
                "polynomial": str(value),
                "total_degree": int(value.total_degree()),
                "degree_r2": int(value.degree(x)),
                "degree_r3": int(value.degree(y)),
                "term_count": len(value.dict()),
                "factorization": implementation.factor_record(value, x, y),
            }
            for value in plane_basis
        ],
        "plane_basis_sha256": hashlib.sha256(plane_text.encode()).hexdigest(),
        "limitations": [
            "The lower-recurrence determinant-zero divisor is a separate chart.",
            "The zero-discriminant component is excluded by an exact open scalar.",
            "A plane projection need not be birational or irreducible.",
            "A single finite-field equation is not an equation over Q.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    })
    result["record_sha256"] = implementation.canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        result = compute(arguments.prime, arguments.census, arguments.output)
        print(json.dumps({
            key: result.get(key)
            for key in (
                "status", "prime", "saturated_krull_dimension",
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
