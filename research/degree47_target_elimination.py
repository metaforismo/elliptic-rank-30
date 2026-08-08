#!/usr/bin/env python3
"""Exact elimination for the generic maximal degree-(4,7) target system.

Under the squarefree/coprime Mason equality case, normalize Q(0)=1 and use

    Q^2-v^2L^3=v^3-Sv^2+3v+1.

The Wronskian identity recursively determines Q from L.  After setting

    L=A(1 + b v^2 + c v^3 + d v^4),

(the linear coefficient is forced to zero), the v^3 target coefficient gives

    c=2/11-3b.

Three exact equations in b,d remain.  This script computes their Groebner
elimination over Q, enumerates every rational pair, applies the leading
cube-class condition, and reconstructs/verifies every full rational identity.

SymPy is used only with exact QQ polynomial arithmetic.  The generated JSON
stores the full basis, factors, candidates, and exact reconstructions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp


b, d = sp.symbols("b d")
c = sp.Rational(2, 11) - 3 * b


def reduced_equations() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
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
    return tuple(sp.factor(expression) for expression in (e8, e9, e10))


def q_coefficients(b_value: sp.Expr, d_value: sp.Expr) -> list[sp.Expr]:
    c_value = sp.Rational(2, 11) - 3 * b_value
    q0 = sp.Integer(1)
    q1 = sp.Rational(3, 2)
    q2 = 4 * b_value
    q3 = (9 * b_value + 11 * c_value) / 4
    q4 = (32 * b_value**2 + 27 * c_value + 28 * d_value) / 12
    q5 = (9 * b_value**2 + 67 * b_value * c_value + 36 * d_value) / 16
    q6 = (
        9 * b_value * c_value
        + 11 * c_value**2
        + 32 * b_value * d_value
    ) / 8
    q7 = (
        -9 * b_value**3
        - 3 * b_value**2 * c_value
        + 108 * b_value * d_value
        + 54 * c_value**2
        + 232 * c_value * d_value
    ) / 96
    return [sp.factor(value) for value in (q0, q1, q2, q3, q4, q5, q6, q7)]


def primitive_integer_poly(expression: sp.Expr, *variables: sp.Symbol) -> sp.Poly:
    polynomial = sp.Poly(sp.together(expression), *variables, domain=sp.QQ)
    _, primitive = polynomial.clear_denoms(convert=True)
    _, primitive = primitive.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def rational_cube_root(value: sp.Rational) -> sp.Rational | None:
    value = sp.Rational(value)
    if value == 0:
        return sp.Integer(0)

    def integer_cube_root(number: int) -> int | None:
        sign = -1 if number < 0 else 1
        target = abs(number)
        lo, hi = 0, 1
        while hi**3 < target:
            hi *= 2
        while lo <= hi:
            mid = (lo + hi) // 2
            cube = mid**3
            if cube == target:
                return sign * mid
            if cube < target:
                lo = mid + 1
            else:
                hi = mid - 1
        return None

    numerator = integer_cube_root(int(value.p))
    denominator = integer_cube_root(int(value.q))
    if numerator is None or denominator is None:
        return None
    return sp.Rational(numerator, denominator)


def rational_roots(polynomial: sp.Poly) -> list[sp.Rational]:
    roots = sp.ground_roots(polynomial)
    return sorted(sp.Rational(root) for root in roots)


def common_d_roots(equations: Iterable[sp.Expr], b_value: sp.Rational) -> list[sp.Rational]:
    specialized: list[sp.Poly] = []
    for equation in equations:
        polynomial = primitive_integer_poly(equation.subs(b, b_value), d)
        if polynomial.is_zero:
            continue
        specialized.append(polynomial)
    if not specialized:
        raise AssertionError("all specialized equations vanished")
    gcd_polynomial = specialized[0]
    for polynomial in specialized[1:]:
        gcd_polynomial = sp.gcd(gcd_polynomial, polynomial)
    return rational_roots(gcd_polynomial)


def rational_string(value: sp.Expr) -> str:
    value = sp.factor(value)
    if value.is_Rational:
        rational = sp.Rational(value)
        return str(rational.p) if rational.q == 1 else f"{rational.p}/{rational.q}"
    return str(value)


def verify_full_identity(
    A: sp.Rational,
    b_value: sp.Rational,
    d_value: sp.Rational,
) -> dict:
    v = sp.symbols("v")
    c_value = sp.Rational(2, 11) - 3 * b_value
    qs = q_coefficients(b_value, d_value)
    L = sp.expand(A * (1 + b_value * v**2 + c_value * v**3 + d_value * v**4))
    Q = sp.expand(sum(qs[index] * v**index for index in range(8)))
    residual = sp.Poly(sp.expand(Q**2 - v**2 * L**3), v, domain=sp.QQ)
    if residual.degree() != 3:
        raise AssertionError("reconstructed residual is not cubic")
    coefficients = {degree: residual.nth(degree) for degree in range(4)}
    if coefficients[0] != 1 or coefficients[1] != 3 or coefficients[3] != 1:
        raise AssertionError("reconstructed residual has wrong fixed coefficients")
    S = -coefficients[2]
    expected = sp.Poly(v**3 - S * v**2 + 3 * v + 1, v, domain=sp.QQ)
    if residual != expected:
        raise AssertionError("full target identity failed")
    return {
        "A": rational_string(A),
        "L": rational_string(L),
        "Q": rational_string(Q),
        "S": rational_string(S),
        "b": rational_string(b_value),
        "c": rational_string(c_value),
        "d": rational_string(d_value),
    }


def compute_elimination() -> dict:
    equations = reduced_equations()
    integer_equations = [primitive_integer_poly(expression, d, b) for expression in equations]
    basis = sp.groebner(integer_equations, d, b, order="lex", domain=sp.QQ)
    basis_polynomials = [primitive_integer_poly(poly.as_expr(), d, b) for poly in basis.polys]
    univariate = [poly for poly in basis_polynomials if not poly.as_expr().has(d)]
    if not univariate:
        raise AssertionError("Groebner basis has no b-elimination polynomial")
    elimination = min(univariate, key=lambda poly: poly.degree())
    factor_constant, factors = sp.factor_list(elimination.as_expr())

    b_roots = rational_roots(elimination)
    pairs: list[tuple[sp.Rational, sp.Rational]] = []
    for b_value in b_roots:
        for d_value in common_d_roots(equations, b_value):
            if all(sp.factor(expression.subs({b: b_value, d: d_value})) == 0 for expression in equations):
                pairs.append((b_value, d_value))
    pairs = sorted(set(pairs))

    leading_candidates: list[dict] = []
    full_solutions: list[dict] = []
    for b_value, d_value in pairs:
        c_value = sp.Rational(2, 11) - 3 * b_value
        q7 = sp.Rational(q_coefficients(b_value, d_value)[7])
        record = {
            "b": rational_string(b_value),
            "c": rational_string(c_value),
            "d": rational_string(d_value),
            "q7": rational_string(q7),
        }
        if d_value == 0 or q7 == 0:
            record["leading_status"] = "degenerate_degree"
            leading_candidates.append(record)
            continue
        A_cubed = sp.factor(q7**2 / d_value**3)
        record["A_cubed"] = rational_string(A_cubed)
        A = rational_cube_root(sp.Rational(A_cubed))
        if A is None or A == 0:
            record["leading_status"] = "not_a_rational_cube"
            leading_candidates.append(record)
            continue
        record["leading_status"] = "compatible"
        record["A"] = rational_string(A)
        leading_candidates.append(record)
        full_solutions.append(verify_full_identity(A, b_value, d_value))

    result = {
        "certificate_id": "degree47-target-elimination-v1",
        "claim_status": "proved exact symbolic computation",
        "elimination_polynomial": str(elimination.as_expr()),
        "elimination_polynomial_degree": elimination.degree(),
        "factor_constant": str(factor_constant),
        "factorization": [
            {"factor": str(sp.factor(factor)), "multiplicity": multiplicity}
            for factor, multiplicity in factors
        ],
        "full_rational_solutions": full_solutions,
        "groebner_basis": [str(poly.as_expr()) for poly in basis_polynomials],
        "leading_candidates": leading_candidates,
        "normalized_system": {
            "c": "2/11-3b",
            "equations": [str(sp.factor(expression)) for expression in equations],
            "q1": "3/2",
            "wronskian": "2LQ+3vL'Q-2vLQ'=constant",
        },
        "rational_b_roots": [rational_string(root) for root in b_roots],
        "rational_pair_count": len(pairs),
        "result": "rational_solution_found" if full_solutions else "no_rational_solution",
        "schema_version": 1,
        "scope": "generic squarefree/coprime degree-(4,7) target system; repeated-root S values are covered by a separate certificate",
        "software": {
            "engine": "SymPy exact QQ Groebner basis and polynomial factorization",
            "sympy_version": sp.__version__,
        },
        "theorem": "Every generic rational degree-(4,7) target identity appears in full_rational_solutions; an empty list proves nonexistence in this branch.",
    }
    unhashed = json.dumps(result, sort_keys=True, separators=(",", ":")).encode("utf-8")
    result["record_sha256"] = hashlib.sha256(unhashed).hexdigest()
    return result


def canonical_payload(result: dict) -> str:
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--expect-result", choices=("no_rational_solution", "rational_solution_found"))
    args = parser.parse_args()
    result = compute_elimination()
    if args.expect_result and result["result"] != args.expect_result:
        raise SystemExit(f"expected {args.expect_result}, obtained {result['result']}")
    payload = canonical_payload(result)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check and args.check.read_text(encoding="utf-8") != payload:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print(
        "VERIFIED degree-(4,7) target elimination",
        f"result={result['result']}",
        f"pairs={result['rational_pair_count']}",
        f"full={len(result['full_rational_solutions'])}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
