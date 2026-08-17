#!/usr/bin/env python3
"""Exact endpoint-reduced section incidence for one split IV surface over F_p.

For a fixed normalized I12+I4+IV surface, write

    c4=(t-1)^2*A,  c6=(t-1)^2*B,
    X=(t-1)*W,     Y=t*(t-1)*Z,
    W=e0*P*D^2+t*U,

where D=d0+d1*t+t^2 and P is the certified I4 square jet.  The section identity
then becomes Z^2=H with deg(H)=20.  The constant and leading square equations
are built into the split I4 and smooth-infinity parametrizations.  Since the
constant coefficient z0 is saturated as nonzero, z1,...,z9 are eliminated
recursively.  The remaining coefficient equations, the chosen split-IV branch,
and one Rabinowitsch saturation define an ideal in nine geometric variables.

The script computes its exact finite-field Groebner basis.  If the quotient is
zero-dimensional and reasonably small, it also constructs multiplication
matrices and tests deterministic primitive coordinates.  It makes no lift or
rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import time
import traceback
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def polynomial_record(polynomial) -> dict[str, object]:
    return {
        "degree": int(polynomial.total_degree()) if polynomial else None,
        "term_count": len(polynomial.dict()) if polynomial else 0,
        "sha256": hashlib.sha256(str(polynomial).encode()).hexdigest(),
    }


def factor_record(polynomial) -> list[dict[str, object]]:
    return [
        {
            "factor": str(factor),
            "degree": int(factor.degree()),
            "exponent": int(exponent),
        }
        for factor, exponent in polynomial.factor()
    ]


def monomial_from_exponents(ring, exponents: tuple[int, ...]):
    value = ring.one()
    for generator, exponent in zip(ring.gens(), exponents, strict=True):
        if exponent:
            value *= generator**exponent
    return value


def divisible_by_any(
    exponents: tuple[int, ...],
    leading_exponents: list[tuple[int, ...]],
) -> bool:
    return any(
        all(left >= right for left, right in zip(exponents, lead, strict=True))
        for lead in leading_exponents
    )


def standard_monomials(
    ring,
    leading_exponents: list[tuple[int, ...]],
    dimension: int,
) -> list[tuple[int, ...]]:
    variable_count = ring.ngens()
    zero = (0,) * variable_count
    heap: list[tuple[int, tuple[int, ...]]] = [(0, zero)]
    queued = {zero}
    result: list[tuple[int, ...]] = []
    while heap and len(result) < dimension:
        _degree, exponents = heapq.heappop(heap)
        if divisible_by_any(exponents, leading_exponents):
            continue
        result.append(exponents)
        for index in range(variable_count):
            neighbour = list(exponents)
            neighbour[index] += 1
            neighbour_tuple = tuple(neighbour)
            if neighbour_tuple in queued:
                continue
            queued.add(neighbour_tuple)
            heapq.heappush(
                heap,
                (sum(neighbour_tuple), neighbour_tuple),
            )
    if len(result) != dimension:
        raise AssertionError((
            "failed to enumerate the declared quotient basis",
            len(result),
            dimension,
        ))
    return result


def multiplication_matrix(
    ring,
    ideal,
    basis_exponents: list[tuple[int, ...]],
    multiplier,
):
    field = ring.base_ring()
    index = {exponents: position for position, exponents in enumerate(basis_exponents)}
    size = len(basis_exponents)
    matrix = Matrix(field, size, size)
    for column, exponents in enumerate(basis_exponents):
        monomial = monomial_from_exponents(ring, exponents)
        normal = ideal.reduce(multiplier * monomial)
        for normal_exponents, coefficient in normal.dict().items():
            normal_tuple = tuple(int(value) for value in normal_exponents)
            if normal_tuple not in index:
                raise AssertionError((
                    "normal form contains a nonstandard monomial",
                    normal_tuple,
                ))
            matrix[index[normal_tuple], column] = coefficient
    return matrix


def load_candidate(path: Path, index: int) -> tuple[dict[str, object], dict[str, object]]:
    certificate = json.loads(path.read_text(encoding="utf-8"))
    candidates = certificate["candidates"]
    if index < 0 or index >= len(candidates):
        raise IndexError((index, len(candidates)))
    return certificate, candidates[index]


def compute(
    *,
    certificate_path: Path,
    candidate_index: int,
    i4_root: int,
    iv_root: int,
    multiplication_cap: int,
    output: Path,
) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    certificate, candidate = load_candidate(certificate_path, candidate_index)
    prime = int(certificate["prime"])
    field = GF(prime)
    parameters = candidate["parameters"]
    e0 = field(parameters["e0"])
    p0 = field(parameters["p0"])
    p1 = field(parameters["p1"])
    p2 = field(parameters["p2"])
    p3 = field(parameters["p3"])
    rho = field(i4_root)
    eta = field(iv_root)

    valid_i4_roots = {
        int(value) for value in candidate["split_tangent_checks"]["i4_square_roots"]
    }
    valid_iv_roots = {
        int(value) for value in candidate["split_tangent_checks"]["iv_square_roots"]
    }
    if int(rho) not in valid_i4_roots:
        raise ValueError(("invalid I4 root", int(rho), sorted(valid_i4_roots)))
    if int(eta) not in valid_iv_roots:
        raise ValueError(("invalid IV root", int(eta), sorted(valid_iv_roots)))
    if rho**2 != -3 * e0 * p0:
        raise AssertionError("I4 tangent root failed its square equation")
    b_at_one = field(
        candidate["exact_fibre_checks"]["c6_quadratic_coefficient_at_one"]
    )
    if eta**2 != -2 * b_at_one:
        raise AssertionError("IV tangent root failed its square equation")

    base_names = ("d0", "d1", "q", "u0", "u1", "u2", "u3", "u4", "u5")
    base = PolynomialRing(field, names=base_names, order="degrevlex")
    d0, d1, q, u0, u1, u2, u3, u4, u5 = base.gens()
    fraction_field = base.fraction_field()
    time_ring = PolynomialRing(base, names=("t",))
    t = time_ring.gen()

    c4 = time_ring([
        field(value) for value in candidate["c4_coefficients_ascending"]
    ])
    c6 = time_ring([
        field(value) for value in candidate["c6_coefficients_ascending"]
    ])
    divisor = (t - 1) ** 2
    A, remainder_a = c4.quo_rem(divisor)
    B, remainder_b = c6.quo_rem(divisor)
    if remainder_a or remainder_b:
        raise AssertionError("the candidate is not divisible by (t-1)^2")

    D = d0 + d1 * t + t**2
    P = p0 + p1 * t + p2 * t**2 + p3 * t**3
    u6 = q**2 + 2 - e0 * p3
    U = u0 + u1*t + u2*t**2 + u3*t**3 + u4*t**4 + u5*t**5 + u6*t**6
    W = e0 * P * D**2 + t * U
    X = (t - 1) * W
    numerator = (t - 1) * (W**3 - 3 * A * W * D**4) - 2 * B * D**6
    H, remainder_h = numerator.quo_rem(t**2)
    if remainder_h:
        raise AssertionError("the I4-adapted numerator is not divisible by t^2")
    if H.degree() > 20:
        raise AssertionError(("unexpected H degree", H.degree()))
    coefficients = [base(H[index]) for index in range(21)]

    node_factor = e0 * p0 * d0**2 - u0
    z0 = rho * d0 * node_factor
    z10 = q * (q**2 + 3)
    if coefficients[0] != z0**2:
        raise AssertionError(("constant square equation failed", coefficients[0] - z0**2))
    if coefficients[20] != z10**2:
        raise AssertionError(("leading square equation failed", coefficients[20] - z10**2))

    z = [fraction_field(z0)]
    two_z0 = fraction_field(2 * z0)
    for degree in range(1, 10):
        correction = sum(
            (z[index] * z[degree - index] for index in range(1, degree)),
            fraction_field.zero(),
        )
        z.append((fraction_field(coefficients[degree]) - correction) / two_z0)
    z.append(fraction_field(z10))

    rational_residuals = []
    residual_labels = []
    for degree in range(10, 20):
        square_coefficient = sum(
            (
                z[index] * z[degree - index]
                for index in range(max(0, degree - 10), min(10, degree) + 1)
            ),
            fraction_field.zero(),
        )
        residual = square_coefficient - fraction_field(coefficients[degree])
        if residual:
            rational_residuals.append(residual)
            residual_labels.append(f"square_coefficient_{degree}")

    leading_residual = z10**2 - fraction_field(coefficients[20])
    if leading_residual:
        raise AssertionError("leading residual reappeared after recursion")

    D1 = d0 + d1 + 1
    iv_residual = sum(z, fraction_field.zero()) - fraction_field(eta * D1**3)
    if iv_residual:
        rational_residuals.append(iv_residual)
        residual_labels.append("chosen_IV_tangent")

    residual_numerators = []
    residual_denominators = []
    for residual in rational_residuals:
        numerator_polynomial = base(residual.numerator())
        denominator_polynomial = base(residual.denominator())
        if numerator_polynomial:
            residual_numerators.append(numerator_polynomial)
            residual_denominators.append(denominator_polynomial)

    resultant_d_x = base(D.resultant(X))
    saturation = base(
        d0 * D1 * node_factor * (q**2 + 3) * resultant_d_x
    )
    if not saturation:
        raise AssertionError("the saturation product vanished identically")

    ring_names = base_names + ("inv",)
    ring = PolynomialRing(field, names=ring_names, order="degrevlex")
    ring_generators = ring.gens()
    embedding = base.hom(ring_generators[: len(base_names)], ring)
    inv = ring_generators[-1]
    equations = [embedding(value) for value in residual_numerators]
    equations.append(inv * embedding(saturation) - 1)

    initial = {
        "status": "groebner_started",
        "truth_status": (
            f"exact endpoint-reduced section incidence over F_{prime}; "
            "no p-adic, characteristic-zero, or rank-30 conclusion"
        ),
        "sage_version": str(sage_version),
        "prime": prime,
        "surface_certificate_sha256": certificate["certificate_sha256"],
        "candidate_index": candidate_index,
        "candidate_record_sha256": candidate["record_sha256"],
        "i4_root": int(rho),
        "iv_root": int(eta),
        "base_variables": list(base_names),
        "ring_variables": list(ring_names),
        "residual_labels": residual_labels,
        "residual_equations": [polynomial_record(value) for value in residual_numerators],
        "residual_denominators": [polynomial_record(value) for value in residual_denominators],
        "saturation": polynomial_record(saturation),
        "resultant_D_X": polynomial_record(resultant_d_x),
        "multiplication_cap": multiplication_cap,
    }
    output.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")

    started = time.time()
    ideal = ring.ideal(equations)
    std = singular_function("std")
    dim_function = singular_function("dim")
    vdim_function = singular_function("vdim")
    singular_basis = std(ideal)
    elapsed = time.time() - started
    dimension = int(dim_function(singular_basis))
    quotient_dimension = int(vdim_function(singular_basis)) if dimension == 0 else None
    basis = [ring(polynomial) for polynomial in singular_basis]
    basis_text = "\n".join(str(polynomial) for polynomial in basis) + "\n"
    (output.parent / "groebner-basis.txt").write_text(basis_text)

    result = dict(initial)
    result.update({
        "status": "completed",
        "elapsed_seconds": elapsed,
        "equation_count": len(equations),
        "equation_degrees": [int(value.total_degree()) for value in equations],
        "groebner_basis_count": len(basis),
        "groebner_basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        "krull_dimension": dimension,
        "quotient_dimension": quotient_dimension,
        "multiplication_coordinates": [],
    })

    if dimension == 0 and quotient_dimension is not None and quotient_dimension <= multiplication_cap:
        leading_exponents = [
            tuple(int(value) for value in polynomial.lm().exponents()[0])
            for polynomial in basis
            if polynomial
        ]
        basis_exponents = standard_monomials(
            ring,
            leading_exponents,
            quotient_dimension,
        )
        groebner_ideal = ring.ideal(basis)
        d0_s, d1_s, q_s, u0_s, u1_s, u2_s, u3_s, u4_s, u5_s, _inv_s = ring.gens()
        coordinate_candidates = [
            ("q_squared", q_s**2),
            ("q", q_s),
            ("d0", d0_s),
            ("d1", d1_s),
            ("u0", u0_s),
            (
                "deterministic_linear_form",
                q_s + 2*d0_s + 3*d1_s + 4*u0_s + 5*u1_s
                + 6*u2_s + 7*u3_s + 8*u4_s + 9*u5_s,
            ),
        ]
        for name, coordinate in coordinate_candidates:
            matrix = multiplication_matrix(
                ring,
                groebner_ideal,
                basis_exponents,
                coordinate,
            )
            characteristic = matrix.charpoly("T")
            try:
                minimal = matrix.minimal_polynomial()
            except AttributeError:
                minimal = matrix.minpoly()
            squarefree = minimal.gcd(minimal.derivative()).degree() == 0
            matrix_raw = json.dumps(
                [[int(value) for value in row] for row in matrix.rows()],
                separators=(",", ":"),
            ).encode()
            result["multiplication_coordinates"].append({
                "name": name,
                "expression": str(coordinate),
                "matrix_sha256": hashlib.sha256(matrix_raw).hexdigest(),
                "minimal_polynomial": str(minimal),
                "minimal_polynomial_degree": int(minimal.degree()),
                "minimal_polynomial_factors": factor_record(minimal),
                "characteristic_polynomial": str(characteristic),
                "characteristic_polynomial_factors": factor_record(characteristic),
                "separates_quotient": int(minimal.degree()) == quotient_dimension,
                "minimal_polynomial_squarefree": bool(squarefree),
            })
            if int(minimal.degree()) == quotient_dimension and squarefree:
                result["primitive_etale_coordinate"] = name
                break
        result["standard_monomial_basis_sha256"] = canonical_hash(
            [list(exponents) for exponents in basis_exponents]
        )
    elif dimension == 0:
        result["multiplication_skip_reason"] = (
            f"quotient dimension {quotient_dimension} exceeds cap {multiplication_cap}"
        )

    result["record_sha256"] = canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--candidate-index", type=int, required=True)
    parser.add_argument("--i4-root", type=int, required=True)
    parser.add_argument("--iv-root", type=int, required=True)
    parser.add_argument("--multiplication-cap", type=int, default=384)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    try:
        result = compute(
            certificate_path=arguments.certificate,
            candidate_index=arguments.candidate_index,
            i4_root=arguments.i4_root,
            iv_root=arguments.iv_root,
            multiplication_cap=arguments.multiplication_cap,
            output=arguments.output,
        )
        print(json.dumps({
            key: result.get(key)
            for key in (
                "status",
                "candidate_index",
                "i4_root",
                "iv_root",
                "krull_dimension",
                "quotient_dimension",
                "primitive_etale_coordinate",
                "elapsed_seconds",
                "record_sha256",
            )
        }, indent=2, sort_keys=True), flush=True)
    except Exception as exc:
        error = {
            "status": "error",
            "truth_status": (
                "implementation or computation failure; no mathematical conclusion"
            ),
            "error": repr(exc),
            "traceback": traceback.format_exc(),
        }
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(error, indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(error, indent=2, sort_keys=True), flush=True)
        raise


if __name__ == "__main__":
    main()
