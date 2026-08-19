#!/usr/bin/env python3
"""Exact rational obstruction on the determinant-zero additive-IV chart.

The logarithmic-derivative reduction solves a 3x3 affine-linear system in
(m0,r1,r0).  Its determinant is

    D = -2*k^2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3^2 + 148*r3 - 320.

The regular chart inverts D.  This script treats D=0 directly, without clearing
that determinant and without inferring anything from the regular formulas.

After imposing D=0, two lower equations express r1 and r0 in terms of the free
parameter u=m0.  The third equation is exactly the plane quintic G5(k,r3)=0.
The next differential coefficient is affine-linear in u, a8*u+b8.

* On a8=0, consistency also requires b8=0.  Exact resultants with G5 have
  coprime eliminants in r3, proving this affine branch empty over Qbar.
* On a8!=0, substitute u=-b8/a8.  The next two exact coefficients give plane
  equations H7,H6.  The gcd of Res_k(G5,H7) and Res_k(G5,H6) is an irreducible
  degree-9 polynomial in r3, hence it has no rational root.

Therefore the affine determinant-zero chart contains no Q-rational solution of
the complete differential system.  This does not cover projective coordinate
boundaries outside the declared affine chart and does not change the rank-29
lower bound.
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


def primitive(poly, *variables):
    value = sp.Poly(poly, *variables, domain=sp.QQ)
    _, value = value.clear_denoms(convert=True)
    _, value = value.primitive()
    if value.LC() < 0:
        value = -value
    return value


def primitive_resultant(first, second, eliminated, remaining):
    value = sp.resultant(first, second, eliminated)
    numerator, _denominator = sp.cancel(value).as_numer_denom()
    return primitive(numerator, remaining).monic()


def build() -> dict[str, object]:
    t = sp.symbols("t")
    r0, r1, r2, r3, k = sp.symbols("r0 r1 r2 r3 k")
    m0, m1, m2, u = sp.symbols("m0 m1 m2 u")

    R = r0 + r1*t + r2*t**2 + r3*t**3 + t**4
    P = t*R
    Q = sp.diff(P, t) + 3*R
    K = k - 12*t
    M = m0 + m1*t + m2*t**2 + 80*t**3
    N = sp.expand(Q**2 + M*P)
    T = sp.expand(
        N*K*((t-1)*Q - 2*P)
        - 3*P*(t-1)*(sp.diff(N, t)*K - 2*N*sp.diff(K, t))
    )
    E = sp.expand(
        Q*T*K - 2*P*sp.diff(T, t)*K
        + 8*P*T*sp.diff(K, t) - (t-1)*N**2*K**2
    )
    quotient, remainder = sp.div(sp.Poly(E, t), sp.Poly(P, t))
    if not remainder.is_zero:
        raise AssertionError("E lost its forced P divisor")
    qbar = sp.Poly(quotient.as_expr(), t)
    coefficient = {
        degree: qbar.coeff_monomial(t**degree)
        for degree in range(6, 14)
    }

    m2_formula = 8*(-k + 2*r3 + 4)
    m1_formula = -(
        k**2 + 14*k*r3 + 4*k + 24*r2
        + 13*r3**2 - 56*r3 - 80
    )/3
    top_substitution = {m2: m2_formula, m1: m1_formula}
    lower_equations = [
        sp.expand(coefficient[degree].subs(top_substitution))
        for degree in (11, 10, 9)
    ]
    matrix, rhs = sp.linear_eq_to_matrix(
        lower_equations, [m0, r1, r0]
    )
    D = (
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    determinant_ratio = sp.factor(matrix.det()/D)
    if not determinant_ratio.is_Rational or determinant_ratio == 0:
        raise AssertionError("unexpected lower-system determinant")

    r2_on_D = sp.solve(sp.Eq(D, 0), r2)[0]
    equations_on_D = [
        sp.cancel(value.subs(r2, r2_on_D))
        for value in lower_equations
    ]
    singular_solution = sp.solve(
        [
            equations_on_D[0].subs(m0, u),
            equations_on_D[1].subs(m0, u),
        ],
        [r1, r0],
        dict=True,
        simplify=False,
    )
    if len(singular_solution) != 1:
        raise AssertionError("the singular lower chart was not one-dimensional")
    singular_solution = singular_solution[0]

    compatibility = sp.cancel(
        equations_on_D[2].subs({
            m0: u,
            r1: singular_solution[r1],
            r0: singular_solution[r0],
        })
    )
    compatibility_numerator, compatibility_denominator = (
        compatibility.as_numer_denom()
    )
    if compatibility_numerator.has(u):
        raise AssertionError("the singular compatibility still depends on u")
    G5 = primitive(compatibility_numerator, k, r3).monic()
    compatibility_ratio = sp.factor(
        compatibility_numerator/G5.as_expr()
    )
    if compatibility_ratio.has(k, r3, u):
        raise AssertionError("compatibility is not a scalar multiple of G5")

    singular_substitution = {
        m2: m2_formula,
        m1: m1_formula,
        r2: r2_on_D,
        m0: u,
        r1: singular_solution[r1],
        r0: singular_solution[r0],
    }

    def reduce_mod_G5(degree: int):
        value = sp.cancel(
            coefficient[degree].subs(singular_substitution)
        )
        numerator, _denominator = value.as_numer_denom()
        numerator = primitive(numerator, k, r3, u)
        coefficient_field = sp.QQ.frac_field(r3, u)
        remainder = sp.rem(
            sp.Poly(numerator.as_expr(), k, domain=coefficient_field),
            sp.Poly(G5.as_expr(), k, domain=coefficient_field),
        )
        remainder_numerator, _remainder_denominator = (
            sp.cancel(remainder.as_expr()).as_numer_denom()
        )
        return primitive(remainder_numerator, k, r3, u)

    residual8 = reduce_mod_G5(8)
    residual7 = reduce_mod_G5(7)
    residual6 = reduce_mod_G5(6)
    if residual8.degree(u) != 1:
        raise AssertionError("the degree-8 singular equation is not affine in u")
    a8 = sp.Poly(residual8.as_expr(), u).coeff_monomial(u)
    b8 = sp.Poly(residual8.as_expr(), u).coeff_monomial(1)
    if sp.gcd(primitive(a8, k, r3), G5).degree() != 0:
        raise AssertionError("a8 vanishes identically on a G5 component")

    # Branch a8=0: residual8 also requires b8=0.
    resultant_a = primitive_resultant(
        G5.as_expr(), a8, k, r3
    )
    resultant_b = primitive_resultant(
        G5.as_expr(), b8, k, r3
    )
    zero_branch_gcd = sp.gcd(resultant_a, resultant_b).monic()
    zero_branch_empty = zero_branch_gcd.degree() == 0
    if not zero_branch_empty:
        raise AssertionError(
            "the a8=0 branch has unresolved affine points"
        )

    # Branch a8!=0: solve for u and use the next two coefficients.
    u_formula = sp.cancel(-b8/a8)

    def condition_after_u(residual):
        value = sp.cancel(residual.as_expr().subs(u, u_formula))
        numerator, _denominator = value.as_numer_denom()
        coefficient_field = sp.QQ.frac_field(r3)
        remainder = sp.rem(
            sp.Poly(numerator, k, domain=coefficient_field),
            sp.Poly(G5.as_expr(), k, domain=coefficient_field),
        )
        remainder_numerator, _remainder_denominator = (
            sp.cancel(remainder.as_expr()).as_numer_denom()
        )
        return primitive(remainder_numerator, k, r3).monic()

    H7 = condition_after_u(residual7)
    H6 = condition_after_u(residual6)
    resultant7 = primitive_resultant(
        G5.as_expr(), H7.as_expr(), k, r3
    )
    resultant6 = primitive_resultant(
        G5.as_expr(), H6.as_expr(), k, r3
    )
    H9 = sp.gcd(resultant7, resultant6).monic()
    factorization = sp.factor_list(H9.as_expr())
    factor_degrees = [
        {"degree": int(sp.degree(factor, r3)), "exponent": int(exponent)}
        for factor, exponent in factorization[1]
    ]
    has_linear_factor = any(item["degree"] == 1 for item in factor_degrees)
    if H9.degree() != 9 or has_linear_factor:
        raise AssertionError("the nonzero-a8 eliminant does not exclude Q-points")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_singular_chart_rational_obstruction",
        "truth_status": (
            "CERTIFIED absence of Q-rational solutions on the declared affine "
            "D=0 logarithmic-derivative chart; projective coordinate boundaries, "
            "the semistable branch, and rank 30 remain unresolved"
        ),
        "lower_determinant": str(D),
        "lower_determinant_ratio": str(determinant_ratio),
        "r2_on_determinant": str(sp.factor(r2_on_D)),
        "singular_lower_parametrization": {
            "free_parameter": "u=m0",
            "r1": str(sp.factor(singular_solution[r1])),
            "r0": str(sp.factor(singular_solution[r0])),
        },
        "compatibility_quintic_G5": {
            "polynomial": str(primitive(G5.as_expr(), k, r3).as_expr()),
            "total_degree": int(G5.total_degree()),
            "term_count": len(G5.terms()),
            "compatibility_scalar": str(compatibility_ratio),
            "compatibility_denominator": str(compatibility_denominator),
        },
        "degree8_affine_equation": {
            "a8": str(primitive(a8, k, r3).as_expr()),
            "b8": str(primitive(b8, k, r3).as_expr()),
            "a8_total_degree": int(sp.Poly(a8, k, r3).total_degree()),
            "b8_total_degree": int(sp.Poly(b8, k, r3).total_degree()),
        },
        "a8_zero_branch": {
            "resultant_a_degree": int(resultant_a.degree()),
            "resultant_b_degree": int(resultant_b.degree()),
            "eliminant_gcd": str(zero_branch_gcd.as_expr()),
            "empty_over_algebraic_closure": zero_branch_empty,
        },
        "a8_nonzero_branch": {
            "u_formula": str(sp.factor(u_formula)),
            "H7_total_degree": int(H7.total_degree()),
            "H6_total_degree": int(H6.total_degree()),
            "resultant7_degree": int(resultant7.degree()),
            "resultant6_degree": int(resultant6.degree()),
            "common_eliminant_H9": str(
                primitive(H9.as_expr(), r3).as_expr()
            ),
            "H9_degree": int(H9.degree()),
            "H9_factor_degrees": factor_degrees,
            "H9_has_rational_root": has_linear_factor,
        },
        "conclusion": {
            "a8_zero_branch_empty": zero_branch_empty,
            "a8_nonzero_branch_has_no_Q_point": not has_linear_factor,
            "affine_D_zero_chart_has_no_Q_point": (
                zero_branch_empty and not has_linear_factor
            ),
        },
        "limitations": [
            "The proof concerns the affine logarithmic-derivative coordinates with finite k,r2,r3.",
            "Projective boundary charts must be checked separately before eliminating the whole additive-IV strategy.",
            "No split-square condition or height-79/12 section is needed for this obstruction.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != payload:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
    print(json.dumps({
        "affine_D_zero_chart_has_no_Q_point": payload["conclusion"][
            "affine_D_zero_chart_has_no_Q_point"
        ],
        "H9_degree": payload["a8_nonzero_branch"]["H9_degree"],
        "H9_factor_degrees": payload["a8_nonzero_branch"][
            "H9_factor_degrees"
        ],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
