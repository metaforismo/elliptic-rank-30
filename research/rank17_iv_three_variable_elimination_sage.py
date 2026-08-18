#!/usr/bin/env sage -python
"""Exact modular elimination for the three-variable additive-IV surface model.

The logarithmic-derivative reduction leaves the coordinates (r2,r3,k).  The
remaining lower coefficients are solved over the rational function field, the
residual differential coefficients are converted back to polynomial equations,
and the determinant divisor is saturated by a Rabinowitsch variable.  A lexicographic
Singular basis then eliminates h and k to the (r2,r3) plane.

Every source surface from the supplied finite-field census is independently
mapped into the reduced chart and checked against every residual equation.  A
single modular eliminant is not an equation over Q and this script makes no
section or rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
import traceback
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, vector
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version

HERE = Path(__file__).resolve().parent
BASE_IMPLEMENTATION = HERE / "rank17_iv_log_derivative_three_variable.py"
BAD_PRIMES = {2, 3, 5, 11}


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_base_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_log_derivative_three_variable", BASE_IMPLEMENTATION
    )
    if spec is None or spec.loader is None:
        raise ImportError(BASE_IMPLEMENTATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def construct_residuals(prime: int):
    if prime in BAD_PRIMES:
        raise ValueError(("bad prime for this coordinate chart", prime))
    field = GF(prime)
    base = PolynomialRing(
        field, names=("r2", "r3", "k"), order="degrevlex"
    )
    r2, r3, k = base.gens()
    fraction = base.fraction_field()
    time_ring = PolynomialRing(fraction, "t")
    t = time_ring.gen()
    s = t - 1

    m2 = fraction(8 * (-k + 2*r3 + 4))
    m1 = fraction(-(
        k**2 + 14*k*r3 + 4*k + 24*r2
        + 13*r3**2 - 56*r3 - 80
    )) / 3

    def qbar(*, r0_value=0, r1_value=0, m0_value=0):
        r0f = fraction(r0_value)
        r1f = fraction(r1_value)
        m0f = fraction(m0_value)
        R = r0f + r1f*t + r2*t**2 + r3*t**3 + t**4
        P = t*R
        Q = P.derivative() + 3*R
        K = k - 12*t
        M = m0f + m1*t + m2*t**2 + 80*t**3
        N = Q**2 + M*P
        T = (
            N*K*(s*Q - 2*P)
            - 3*P*s*(N.derivative()*K - 2*N*K.derivative())
        )
        E = (
            Q*T*K - 2*P*T.derivative()*K
            + 8*P*T*K.derivative() - s*N**2*K**2
        )
        quotient, remainder = E.quo_rem(P)
        if remainder:
            raise AssertionError("E lost its forced P factor")
        return quotient

    degrees = (11, 10, 9)
    zero = qbar()
    base_values = vector(fraction, [zero[degree] for degree in degrees])
    columns = []
    for keyword in ("m0_value", "r1_value", "r0_value"):
        arguments = {keyword: 1}
        trial = qbar(**arguments)
        columns.append(vector(
            fraction,
            [trial[degree] - zero[degree] for degree in degrees],
        ))
    matrix = Matrix(fraction, 3, 3, lambda row, column: columns[column][row])
    determinant = matrix.det()
    if not determinant:
        raise AssertionError("lower recurrence determinant vanished identically")
    solution = matrix.solve_right(-base_values)
    m0_value, r1_value, r0_value = solution

    expected_determinant = base(
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    ratio = fraction(determinant) / fraction(expected_determinant)
    if ratio.numerator().degree() != 0 or ratio.denominator().degree() != 0:
        raise AssertionError(("unexpected determinant factor", determinant))

    final = qbar(
        r0_value=r0_value,
        r1_value=r1_value,
        m0_value=m0_value,
    )
    for degree in (13, 12, 11, 10, 9):
        if final[degree]:
            raise AssertionError(("solved coefficient survived", degree, final[degree]))

    residuals = []
    residual_metadata = []
    seen = set()
    for degree in range(8, -1, -1):
        value = final[degree]
        numerator = base(value.numerator())
        denominator = base(value.denominator())
        if not numerator:
            continue
        determinant_power = 0
        while True:
            quotient, remainder = numerator.quo_rem(expected_determinant)
            if remainder:
                break
            numerator = quotient
            determinant_power += 1
        numerator = numerator.monic()
        key = str(numerator)
        if key in seen:
            continue
        seen.add(key)
        residuals.append(numerator)
        residual_metadata.append({
            "differential_degree": degree,
            "total_degree": int(numerator.total_degree()),
            "degree_r2": int(numerator.degree(r2)),
            "degree_r3": int(numerator.degree(r3)),
            "degree_k": int(numerator.degree(k)),
            "term_count": len(numerator.dict()),
            "removed_determinant_power": determinant_power,
            "original_denominator": str(denominator),
            "sha256": hashlib.sha256(str(numerator).encode()).hexdigest(),
        })
    if not residuals:
        raise AssertionError("no residual polynomial survived")
    return {
        "field": field,
        "base": base,
        "residuals": residuals,
        "residual_metadata": residual_metadata,
        "determinant": expected_determinant,
        "determinant_ratio": str(ratio),
        "lower_solution": {
            "m0": str(m0_value),
            "r1": str(r1_value),
            "r0": str(r0_value),
            "m1": str(m1),
            "m2": str(m2),
        },
    }


def replay_source(census_path: Path, construction):
    source = json.loads(census_path.read_text(encoding="utf-8"))
    prime = int(source["prime"])
    module = load_base_module()
    base = construction["base"]
    r2, r3, k = base.gens()
    records = []
    for source_record in source["records"]:
        reduced = module.derive_from_surface(source_record, prime)
        coordinates = reduced["reduced_coordinates"]
        substitution = {
            r2: int(coordinates["r2"]),
            r3: int(coordinates["r3"]),
            k: int(coordinates["k"]),
        }
        values = [int(polynomial(substitution)) for polynomial in construction["residuals"]]
        determinant_value = int(construction["determinant"](substitution))
        if any(values) or determinant_value == 0:
            raise AssertionError((
                "source surface missed the saturated reduced chart",
                source_record["record_sha256"],
                values,
                determinant_value,
            ))
        records.append({
            "source_record_sha256": source_record["record_sha256"],
            "r2": int(coordinates["r2"]),
            "r3": int(coordinates["r3"]),
            "k": int(coordinates["k"]),
            "determinant": determinant_value,
        })
    return {
        "source_record_sha256": source["record_sha256"],
        "source_surface_count": len(source["records"]),
        "all_source_points_satisfy_residuals": True,
        "all_source_points_avoid_determinant": True,
        "points": records,
    }


def compute(prime: int, census_path: Path, output: Path):
    construction = construct_residuals(prime)
    source_replay = replay_source(census_path, construction)
    field = construction["field"]
    base = construction["base"]
    r2_base, r3_base, k_base = base.gens()

    lex = PolynomialRing(
        field, names=("h", "k", "r2", "r3"), order="lex"
    )
    h, k, r2, r3 = lex.gens()
    embedding = base.hom([r2, r3, k], lex)
    equations = [embedding(value) for value in construction["residuals"]]
    equations.append(h*embedding(construction["determinant"]) - 1)

    initial = {
        "schema_version": 1,
        "status": "groebner_started",
        "truth_status": (
            f"exact saturated three-variable additive-IV surface elimination over F_{prime}; "
            "no characteristic-zero, section, or rank-30 conclusion"
        ),
        "prime": prime,
        "sage_version": str(sage_version),
        "variable_order": [str(value) for value in lex.gens()],
        "residual_metadata": construction["residual_metadata"],
        "determinant": str(construction["determinant"]),
        "determinant_ratio": construction["determinant_ratio"],
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
    (output.parent / "three-variable-lex-basis.txt").write_text(basis_text)

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
    (output.parent / "three-variable-plane-basis-r2-r3.txt").write_text(plane_text)

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
                "factorization": factor_record(value, x, y),
            }
            for value in plane_basis
        ],
        "plane_basis_sha256": hashlib.sha256(plane_text.encode()).hexdigest(),
        "limitations": [
            "The determinant-zero divisor is excluded by saturation and remains a separate chart.",
            "A plane projection need not be birational or irreducible.",
            "A single finite-field equation is not an equation over Q.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    })
    result["record_sha256"] = canonical_hash(result)
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
