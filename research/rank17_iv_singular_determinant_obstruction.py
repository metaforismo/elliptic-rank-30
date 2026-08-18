#!/usr/bin/env python3
"""Exact rational obstruction on the determinant-zero IV surface chart.

This script starts from the logarithmic-derivative equation for the normalized
I12+I4+IV surface locus, imposes the singular lower-recurrence determinant D=0,
and proves that no rational point on that chart satisfies all remaining
coefficient equations.

It makes no Mordell-Weil, characteristic-zero section, or rank-30 claim outside
this declared chart.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def primitive_poly(expression, *gens):
    p = sp.Poly(sp.together(expression).as_numer_denom()[0], *gens, domain=sp.QQ)
    _, z = p.clear_denoms(convert=True)
    _, z = z.primitive()
    if z.LC() < 0:
        z = -z
    return z


def coeff_in_u_mod_g(expression, exponent, u, k, field, modulus):
    p = sp.Poly(expression, u, k, domain=field)
    result = sp.Poly(0, k, domain=field)
    for (u_exp, k_exp), coefficient in p.terms():
        if u_exp == exponent:
            result += sp.Poly(coefficient * k**k_exp, k, domain=field)
    return result.rem(modulus)


def build_certificate() -> dict[str, object]:
    t = sp.symbols("t")
    r0, r1, r2, r3, k, m0, m1, m2 = sp.symbols(
        "r0 r1 r2 r3 k m0 m1 m2"
    )
    s = t - 1
    R = r0 + r1*t + r2*t**2 + r3*t**3 + t**4
    P = t*R
    Q = sp.diff(P, t) + 3*R
    K = k - 12*t
    M = m0 + m1*t + m2*t**2 + 80*t**3
    N = sp.expand(Q**2 + M*P)
    T = sp.expand(
        N*K*(s*Q - 2*P)
        - 3*P*s*(sp.diff(N, t)*K - 2*N*sp.diff(K, t))
    )
    E = sp.expand(
        Q*T*K - 2*P*sp.diff(T, t)*K
        + 8*P*T*sp.diff(K, t) - s*N**2*K**2
    )
    qbar, remainder = sp.div(sp.Poly(E, t), sp.Poly(P, t))
    if remainder.as_expr() != 0 or qbar.degree() != 13:
        raise AssertionError("the differential equation lost its exact P divisor")

    m2_formula = 8*(-k + 2*r3 + 4)
    m1_formula = -(
        k**2 + 14*k*r3 + 4*k + 24*r2
        + 13*r3**2 - 56*r3 - 80
    ) / 3
    q = sp.Poly(
        sp.expand(qbar.as_expr().subs({m2: m2_formula, m1: m1_formula})),
        t,
    )
    if q.coeff_monomial(t**13) or q.coeff_monomial(t**12):
        raise AssertionError("the triangular top recurrence failed")

    equations = [sp.expand(q.coeff_monomial(t**d)) for d in (11, 10, 9)]
    matrix, rhs = sp.linear_eq_to_matrix(equations, [m0, r1, r0])
    determinant = sp.factor(matrix.det())
    D = sp.expand(
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    determinant_ratio = sp.cancel(determinant / D)
    if determinant_ratio.free_symbols:
        raise AssertionError("the lower determinant has an unexpected factor")

    r2_on_D = sp.factor(
        (2*k**2 + 37*k*r3 - 52*k + 116*r3**2 - 148*r3 + 320) / 108
    )
    singular_matrix = matrix.subs(r2, r2_on_D)
    singular_rhs = rhs.subs(r2, r2_on_D)
    rank_two_minor = sp.factor(singular_matrix.extract([0, 1], [0, 2]).det())
    if rank_two_minor == 0 or rank_two_minor.free_symbols:
        raise AssertionError("the singular lower system did not have uniform rank two")

    block = singular_matrix.extract([0, 1], [0, 2])
    adjusted_rhs = sp.Matrix([
        singular_rhs[row] - singular_matrix[row, 1]*r1
        for row in (0, 1)
    ])
    m0_solution, r0_solution = map(sp.factor, block.inv()*adjusted_rhs)
    compatibility = sp.factor(
        singular_matrix[2, 0]*m0_solution
        + singular_matrix[2, 1]*r1
        + singular_matrix[2, 2]*r0_solution
        - singular_rhs[2]
    )

    G5 = sp.Poly(
        140*k**5 + 5200*k**4*r3 - 6760*k**4
        + 66380*k**3*r3**2 - 165100*k**3*r3 + 114704*k**3
        + 345650*k**2*r3**3 - 1199040*k**2*r3**2
        + 1669164*k**2*r3 - 777536*k**2
        + 761740*k*r3**4 - 3297100*k*r3**3
        + 7351296*k*r3**2 - 7457920*k*r3 - 1155584*k
        + 592835*r3**5 - 3046960*r3**4 + 9785600*r3**3
        - 10923776*r3**2 + 21043456*r3 + 29003776,
        k, r3, domain=sp.QQ,
    )
    compatibility_ratio = sp.cancel(compatibility / G5.as_expr())
    if compatibility_ratio.free_symbols:
        raise AssertionError("the rank-two compatibility was not exactly G5=0")

    substitutions = {r2: r2_on_D, m0: m0_solution, r0: r0_solution}
    q8_raw = sp.cancel(q.coeff_monomial(t**8).subs(substitutions, simultaneous=True))
    q7_raw = sp.cancel(q.coeff_monomial(t**7).subs(substitutions, simultaneous=True))

    field = sp.QQ.frac_field(r3)
    g = sp.Poly(G5.as_expr(), k, domain=field).monic()
    a = coeff_in_u_mod_g(q8_raw, 1, r1, k, field, g)
    b = coeff_in_u_mod_g(q8_raw, 0, r1, k, field, g)
    if sp.gcd(a, g).degree() != 0:
        raise AssertionError("the q8 slope is not generically invertible on G5")

    q7_coefficients = [
        coeff_in_u_mod_g(q7_raw, exponent, r1, k, field, g)
        for exponent in range(3)
    ]
    obstruction = (
        q7_coefficients[0]*a**2
        - q7_coefficients[1]*b*a
        + q7_coefficients[2]*b**2
    ).rem(g)
    C = primitive_poly(obstruction.as_expr(), k, r3)
    if C.degree(k) != 4 or C.degree(r3) != 13:
        raise AssertionError("unexpected generic q7 obstruction degree")
    if len(sp.factor_list(C.as_expr())[1]) != 1:
        raise AssertionError("the generic q7 obstruction unexpectedly factored")

    resultant = primitive_poly(sp.resultant(G5.as_expr(), C.as_expr(), k), r3)
    resultant_factors = sp.factor_list(resultant.as_expr())
    factor_degrees = [
        [int(sp.Poly(factor, r3).degree()), int(exponent)]
        for factor, exponent in resultant_factors[1]
    ]
    factor_degrees.sort()
    if factor_degrees != [[9, 1], [56, 1]]:
        raise AssertionError(("unexpected generic obstruction factor degrees", factor_degrees))
    if resultant.ground_roots():
        raise AssertionError("the generic singular chart has a rational r3 root")

    A = primitive_poly(a.as_expr(), k, r3)
    B = primitive_poly(b.as_expr(), k, r3)
    resultant_a = primitive_poly(sp.resultant(G5.as_expr(), A.as_expr(), k), r3)
    resultant_b = primitive_poly(sp.resultant(G5.as_expr(), B.as_expr(), k), r3)
    exceptional_gcd = sp.gcd(resultant_a, resultant_b)
    if exceptional_gcd.degree() != 0:
        raise AssertionError("q8 slope and constant have an exceptional common branch")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_singular_determinant_rational_obstruction",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "truth_status": (
            "EXACT rational obstruction for the determinant-zero chart of the "
            "three-variable additive-IV logarithmic-derivative system. It does "
            "not classify the determinant-open chart, construct a section, or "
            "prove or disprove rank 30."
        ),
        "lower_recurrence": {
            "determinant": str(D),
            "determinant_ratio": str(determinant_ratio),
            "rank_two_minor": str(rank_two_minor),
            "r2_on_determinant": str(r2_on_D),
            "m0_in_terms_of_r1": str(m0_solution),
            "r0_in_terms_of_r1": str(r0_solution),
            "compatibility_polynomial": str(G5.as_expr()),
            "compatibility_ratio": str(compatibility_ratio),
        },
        "generic_q8_branch": {
            "q8_slope": str(A.as_expr()),
            "q8_constant": str(B.as_expr()),
            "gcd_q8_slope_G5": "1",
            "q7_obstruction": str(C.as_expr()),
            "q7_obstruction_degrees": {"k": C.degree(k), "r3": C.degree(r3)},
            "resultant_in_r3": str(resultant.as_expr()),
            "resultant_factor_degrees": factor_degrees,
            "resultant_rational_roots": {},
        },
        "exceptional_q8_slope_zero": {
            "resultant_G5_q8_slope": str(resultant_a.as_expr()),
            "resultant_G5_q8_constant": str(resultant_b.as_expr()),
            "resultant_gcd": str(exceptional_gcd.as_expr()),
        },
        "mathematical_consequence": (
            "There is no rational solution of the remaining differential "
            "coefficient equations on D=0: if the q8 slope is nonzero, q7 "
            "forces a resultant with no rational r3 root; if the slope is zero, "
            "q8 would also require its constant term to vanish, but the two "
            "resultants are coprime."
        ),
        "limitations": [
            "Only the determinant-zero chart D=0 is excluded.",
            "The determinant-open rational quintic component is treated by a separate split-square obstruction.",
            "The full additive-IV locus still requires an exhaustive determinant-open decomposition.",
            "No height-79/12 section is constructed.",
            "The unconditional global rank lower bound remains 29.",
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def markdown(payload: dict[str, object]) -> str:
    return "\n".join([
        "# Determinant-zero additive-IV obstruction", "",
        "```text", "SOLVED-A: false", "SOLVED-B: false",
        "unconditional global lower bound: rank E(Q) >= 29", "```", "",
        payload["truth_status"], "", "## Consequence", "",
        payload["mathematical_consequence"], "", "## Exact record", "",
        "```json", json.dumps(payload, indent=2, sort_keys=True), "```", "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(markdown(payload))
    if args.compare:
        committed = json.loads(args.compare.read_text())
        if committed != payload:
            raise AssertionError(f"certificate mismatch: {args.compare}")
    print(json.dumps({
        "certificate_id": payload["certificate_id"],
        "record_sha256": payload["record_sha256"],
        "resultant_factor_degrees": payload["generic_q8_branch"]["resultant_factor_degrees"],
        "resultant_rational_roots": payload["generic_q8_branch"]["resultant_rational_roots"],
        "exceptional_resultant_gcd": payload["exceptional_q8_slope_zero"]["resultant_gcd"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
