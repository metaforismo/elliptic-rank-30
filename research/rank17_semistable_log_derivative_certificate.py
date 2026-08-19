#!/usr/bin/env python3
"""Exact logarithmic-derivative reduction for the semistable rank-17 seed.

This certificate records the sparse structure of the normalized
I12+I4+I3 surface problem.  For

    H=A^3-B^2=kappa*S,
    S=t^4*(t-1)^3*R,
    Q=S',

one has

    Q*H-S*H' = A^2*(Q*A-3*S*A') - B*(Q*B-2*S*B').

If gcd(A,B)=gcd(A,S)=1 on the nondegenerate locus, there is a common
polynomial K such that

    B*K   = Q*A-3*S*A',
    A^2*K = Q*B-2*S*B'.

The endpoint orders force

    Q=t^3*(t-1)^2*Q0,
    K=t^3*(t-1)^2*L,
    M=(A*K^2-Q^2)/S=t^3*(t-1)^2*C,

where deg L=2 with leading coefficient -12 and deg C<=4.  Consequently

    U=Q0^2+t*(t-1)*C*R=A*L^2,

    V=Q0*U*L-3*t*(t-1)*R*(U'*L-2*U*L')=B*L^4,

and the whole surface equation is compressed to

    Q0*V*L-2*t*(t-1)*R*(V'*L-4*V*L')-U^2*L^2=0.

The script proves the generic differential identity and replays every reduced
identity exactly on the committed rational semistable seed.  It does not
construct the missing height-79/12 section or a rank-30 curve.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


C4_COEFFICIENTS = [9, -192, 628, -384, 30, 64, -12, 0, 1]
C6_COEFFICIENTS = [
    -27, 864, -7434, 15488, -14841, 2880, 2316,
    -1152, 99, 96, -18, 0, 1,
]
SEED_PARAMETERS = {
    "e0": 1,
    "e1": 1,
    "p0": -3,
    "p1": 32,
    "p2": 66,
    "p3": 768,
    "q0": -12,
    "q1": -12,
    "q2": 0,
    "s": 0,
}


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def ascending_coefficients(polynomial, variable, length: int | None = None):
    poly = sp.Poly(sp.expand(polynomial), variable, domain=sp.QQ)
    size = poly.degree() + 1 if length is None else length
    return [str(poly.coeff_monomial(variable**index)) for index in range(size)]


def exact_polynomial_quotient(numerator, denominator, variable):
    quotient, remainder = sp.div(
        sp.Poly(sp.expand(numerator), variable, domain=sp.QQ),
        sp.Poly(sp.expand(denominator), variable, domain=sp.QQ),
    )
    if not remainder.is_zero:
        raise AssertionError(("nonzero polynomial remainder", remainder.as_expr()))
    return sp.expand(quotient.as_expr())


def generic_identity() -> dict[str, object]:
    t = sp.symbols("t")
    A = sp.Function("A")(t)
    B = sp.Function("B")(t)
    S = sp.Function("S")(t)
    H = A**3 - B**2
    Q = sp.diff(S, t)
    identity = sp.expand(
        Q*H - S*sp.diff(H, t)
        - A**2*(Q*A - 3*S*sp.diff(A, t))
        + B*(Q*B - 2*S*sp.diff(B, t))
    )
    if sp.simplify(identity) != 0:
        raise AssertionError("generic logarithmic-derivative identity failed")
    return {
        "verified": True,
        "identity": (
            "S'*(A^3-B^2)-S*(A^3-B^2)'="
            "A^2*(S'*A-3*S*A')-B*(S'*B-2*S*B')"
        ),
        "common_syzygy_under_coprimality": [
            "B*K=S'*A-3*S*A'",
            "A^2*K=S'*B-2*S*B'",
        ],
    }


def build() -> dict[str, object]:
    t = sp.symbols("t")
    u = t - 1
    A = sum(sp.Integer(value)*t**index for index, value in enumerate(C4_COEFFICIENTS))
    B = sum(sp.Integer(value)*t**index for index, value in enumerate(C6_COEFFICIENTS))
    R = sp.expand((t + 3)*(t**4 - 6*t**2 + 56*t - 3))
    S = sp.expand(t**4*u**3*R)
    H = sp.factor(A**3 - B**2)
    kappa = sp.cancel(H/S)
    if kappa != 139968:
        raise AssertionError(("unexpected discriminant scalar", kappa))

    if sp.gcd(sp.Poly(A, t), sp.Poly(B, t)).degree() != 0:
        raise AssertionError("A and B are not coprime")
    if sp.gcd(sp.Poly(A, t), sp.Poly(S, t)).degree() != 0:
        raise AssertionError("A and S are not coprime")

    Q = sp.diff(S, t)
    first_numerator = sp.expand(Q*A - 3*S*sp.diff(A, t))
    second_numerator = sp.expand(Q*B - 2*S*sp.diff(B, t))
    K_from_first = exact_polynomial_quotient(first_numerator, B, t)
    K_from_second = exact_polynomial_quotient(second_numerator, A**2, t)
    if sp.expand(K_from_first - K_from_second) != 0:
        raise AssertionError("the two common-syzygy quotients disagree")
    K = sp.factor(K_from_first)

    Q0 = exact_polynomial_quotient(Q, t**3*u**2, t)
    expected_Q0 = sp.expand((7*t - 4)*R + t*u*sp.diff(R, t))
    if sp.expand(Q0 - expected_Q0) != 0:
        raise AssertionError("unexpected endpoint-reduced derivative Q0")
    L = exact_polynomial_quotient(K, t**3*u**2, t)
    if sp.degree(L, t) != 2 or sp.LC(sp.Poly(L, t)) != -12:
        raise AssertionError(("unexpected reduced K polynomial L", L))

    M = exact_polynomial_quotient(A*K**2 - Q**2, S, t)
    C = exact_polynomial_quotient(M, t**3*u**2, t)
    if sp.degree(C, t) > 4:
        raise AssertionError(("unexpected C degree", sp.degree(C, t)))

    U = sp.expand(Q0**2 + t*u*C*R)
    if sp.expand(U - A*L**2) != 0:
        raise AssertionError("U=A*L^2 failed")
    V = sp.expand(
        Q0*U*L
        - 3*t*u*R*(sp.diff(U, t)*L - 2*U*sp.diff(L, t))
    )
    if sp.expand(V - B*L**4) != 0:
        raise AssertionError("V=B*L^4 failed")
    residual = sp.expand(
        Q0*V*L
        - 2*t*u*R*(sp.diff(V, t)*L - 4*V*sp.diff(L, t))
        - U**2*L**2
    )
    if residual != 0:
        raise AssertionError("reduced differential equation did not vanish")

    if sp.expand(B*L - (Q0*A - 3*t*u*R*sp.diff(A, t))) != 0:
        raise AssertionError("first endpoint-reduced syzygy failed")
    if sp.expand(A**2*L - (Q0*B - 2*t*u*R*sp.diff(B, t))) != 0:
        raise AssertionError("second endpoint-reduced syzygy failed")

    e0 = SEED_PARAMETERS["e0"]
    e1 = SEED_PARAMETERS["e1"]
    p0 = SEED_PARAMETERS["p0"]
    q0 = SEED_PARAMETERS["q0"]
    endpoint_l0 = sp.cancel(-4*R.subs(t, 0)/(e0*p0))
    endpoint_l1 = sp.cancel(3*R.subs(t, 1)/(e1*q0))
    if sp.cancel(L.subs(t, 0) - endpoint_l0) != 0:
        raise AssertionError("L(0) endpoint formula failed")
    if sp.cancel(L.subs(t, 1) - endpoint_l1) != 0:
        raise AssertionError("L(1) endpoint formula failed")

    split_i4 = -3*e0*p0
    split_i3 = -3*e1*q0
    if split_i4 != 9 or split_i3 != 36:
        raise AssertionError("unexpected split tangent targets")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_semistable_logarithmic_derivative_reduction",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "truth_status": (
            "CERTIFIED structural logarithmic-derivative reduction and exact replay "
            "on the rational split I12+I4+I3 seed; no height-79/12 section or "
            "rank-30 curve is constructed"
        ),
        "generic_identity": generic_identity(),
        "nondegenerate_seed": {
            "parameters": SEED_PARAMETERS,
            "A_equals_c4_coefficients_ascending": C4_COEFFICIENTS,
            "B_equals_c6_coefficients_ascending": C6_COEFFICIENTS,
            "gcd_A_B": "1",
            "gcd_A_S": "1",
            "discriminant_factorization": str(H),
            "discriminant_scalar_kappa": str(kappa),
            "residual_quintic_R": str(R),
            "R_coefficients_ascending": ascending_coefficients(R, t, 6),
        },
        "endpoint_factorization": {
            "S": "t^4*(t-1)^3*R",
            "Q_equals_S_prime": str(sp.expand(Q)),
            "Q_factorization": "t^3*(t-1)^2*Q0",
            "Q0": str(sp.factor(Q0)),
            "Q0_coefficients_ascending": ascending_coefficients(Q0, t, 7),
            "K": str(K),
            "K_factorization": "t^3*(t-1)^2*L",
            "L": str(sp.factor(L)),
            "L_coefficients_ascending": ascending_coefficients(L, t, 3),
            "L_degree": int(sp.degree(L, t)),
            "L_leading_coefficient": "-12",
            "L_endpoint_formulas": {
                "L(0)=-4*R(0)/(e0*p0)": str(endpoint_l0),
                "L(1)=3*R(1)/(e1*q0)": str(endpoint_l1),
            },
            "M=(A*K^2-Q^2)/S": str(sp.factor(M)),
            "M_factorization": "t^3*(t-1)^2*C",
            "C": str(sp.factor(C)),
            "C_coefficients_ascending": ascending_coefficients(C, t, 5),
            "C_degree": int(sp.degree(C, t)),
            "generic_degree_bound_for_C": 4,
        },
        "compressed_surface_system": {
            "definitions": {
                "Q0": "(7*t-4)*R+t*(t-1)*R'",
                "U": "Q0^2+t*(t-1)*C*R",
                "V": (
                    "Q0*U*L-3*t*(t-1)*R*(U'*L-2*U*L')"
                ),
            },
            "exact_identities_verified": [
                "A*L^2=U",
                "B*L^4=V",
                "B*L=Q0*A-3*t*(t-1)*R*A'",
                "A^2*L=Q0*B-2*t*(t-1)*R*B'",
                (
                    "Q0*V*L-2*t*(t-1)*R*(V'*L-4*V*L')-U^2*L^2=0"
                ),
            ],
            "degrees_on_seed": {
                "R": int(sp.degree(R, t)),
                "Q0": int(sp.degree(Q0, t)),
                "L": int(sp.degree(L, t)),
                "C": int(sp.degree(C, t)),
                "U": int(sp.degree(U, t)),
                "V": int(sp.degree(V, t)),
            },
        },
        "target_section_geometry_not_yet_constructed": {
            "canonical_height": "79/12",
            "intersection_P_dot_O": 2,
            "I12_component_class": 0,
            "I4_component_class": "1 or 3 after reversing the cycle",
            "I3_component_class": "2 or 1 after reversing the cycle",
            "split_tangent_square_targets_on_seed": {
                "I4": split_i4,
                "I3": split_i3,
            },
        },
        "mathematical_consequence": (
            "The six-parameter semistable Hermite--Pade surface problem can be "
            "reformulated using a monic quintic R, a quadratic L with fixed lead "
            "-12, a polynomial C of degree at most 4, two exact divisibilities, "
            "and one endpoint-reduced differential identity."
        ),
        "highest_value_next_action": (
            "Derive the bidirectional coefficient recurrence for (R,L,C), prove "
            "finite-field equivalence with the complete F5/F7 semistable censuses, "
            "then build the section incidence directly in these reduced coordinates."
        ),
        "limitations": [
            "The generic common polynomial K uses the declared coprimality assumptions.",
            "The exact replay proves the reduction on the committed rational seed, not yet scheme-level equivalence for the entire semistable surface locus.",
            "No rational section of height 79/12 is certified.",
            "The unconditional global lower bound remains rank E(Q) >= 29.",
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
        "certificate_id": payload["certificate_id"],
        "C_degree": payload["endpoint_factorization"]["C_degree"],
        "L": payload["endpoint_factorization"]["L"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
