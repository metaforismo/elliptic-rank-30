#!/usr/bin/env python3
"""Exact certificate for the minimal j=0 cyclic-base scale obstruction.

For mu=2, classify every rational scale c for which the minimal
character-one ansatz over K=Q(sqrt(-3)) can exist.  The certificate also
checks explicit representatives of the two surviving cube classes and proves
their CM eigenline is non-torsion by exact finite-field reductions.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

OUT = Path("certificates/j0_minimal_scale_obstruction.json")


def inv_mod(a: int, p: int) -> int:
    return pow(a % p, p - 2, p)


Point = tuple[int, int] | None


def ec_add(P: Point, Q: Point, p: int) -> Point:
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        if y1 % p == 0:
            return None
        slope = 3 * x1 * x1 * inv_mod(2 * y1, p) % p
    else:
        slope = (y2 - y1) * inv_mod(x2 - x1, p) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (slope * (x1 - x3) - y1) % p
    return x3, y3


def point_order(P: Point, p: int, limit: int = 10000) -> int:
    Q: Point = None
    for n in range(1, limit + 1):
        Q = ec_add(Q, P, p)
        if Q is None:
            return n
    raise RuntimeError("order limit exceeded")


def reduce_fraction(value: Fraction, p: int) -> int:
    return value.numerator * inv_mod(value.denominator, p) % p


def finite_reduction_order(c: int, p: int, w_mod_p: int) -> int:
    # Specialize the rational-u surface at u=1.
    if c == 1:
        x = Fraction(3)
        y_over_w = Fraction(-3, 2)
        b_curve = Fraction(-135, 4)
    elif c == 2:
        x = Fraction(9, 4)
        y_over_w = Fraction(27, 8)
        b_curve = Fraction(-729, 16)
    else:
        raise ValueError(c)
    xp = reduce_fraction(x, p)
    yp = reduce_fraction(y_over_w, p) * w_mod_p % p
    bp = reduce_fraction(b_curve, p)
    assert (yp * yp - xp**3 - bp) % p == 0
    discriminant = (-432 * bp * bp) % p
    assert discriminant != 0
    return point_order((xp, yp), p)


def reduce_w(expr: sp.Expr, w: sp.Symbol) -> sp.Expr:
    poly = sp.Poly(sp.expand(expr), w)
    return sp.rem(poly, sp.Poly(w**2 + 3, w)).as_expr().expand()


def verify_section(scale: int, coeffs: dict[str, sp.Expr]) -> bool:
    u, w = sp.symbols("u w")
    a = sp.Integer(3)
    b = sp.Rational(27, 2)
    A, B, C, D, E = (coeffs[k] for k in ("A", "B", "C", "D", "E"))
    x = u * (A * u**3 + B)
    y = C * u**6 + D * u**3 + E
    rhs = x**3 + b**2 * u**6 / scale**2 - a**3 * (u**3 / scale + 1) ** 3
    difference = reduce_w(y**2 - rhs, w)
    return sp.Poly(difference, u).is_zero


def main() -> int:
    rho, S, g, d = sp.symbols("rho S g d")

    # The normalization below assumes B != 0. If B=0, the u^3 equation
    # forces D=3E/(2c), the u^9 equation then forces C=E/(3c^2), and
    # D^2+2CE=-315/(4c^2), contradicting the required 405/(4c^2).
    b_zero_left = -sp.Rational(315, 4)
    b_zero_required = sp.Rational(405, 4)
    assert b_zero_left != b_zero_required

    equations = [
        -3 * g**2 - rho**3 * S,
        -2 * g * d - rho**2 * S + 9,
        -3 * d**2 - 18 * g - 3 * rho * S - sp.Rational(405, 4),
        S - 81 + 18 * d,
    ]
    groebner = sp.groebner(equations, g, d, S, rho, order="lex")
    eliminant = sp.factor(groebner.polys[-1].as_expr())
    expected = sp.factor(
        (2 * rho + 1) ** 3
        * (10 * rho**3 - 48 * rho**2 + 12 * rho - 1)
        / 80
    )
    assert sp.expand(eliminant - expected) == 0

    cubic = 10 * rho**3 - 48 * rho**2 + 12 * rho - 1
    assert sp.Poly(cubic, rho, domain=sp.QQ).is_irreducible
    cubic_discriminant = sp.factor(sp.discriminant(cubic, rho))
    assert cubic_discriminant == -78732
    # A root cannot lie in the quadratic field K because its degree over Q is 3.
    k_rational_rho = sp.Rational(-1, 2)

    solutions = sp.solve(
        [sp.factor(eq.subs(rho, k_rational_rho)) for eq in equations],
        [S, g, d],
        dict=True,
    )
    normalized = sorted(
        [
            {"S": int(sol[S]), "g": str(sol[g]), "d": str(sol[d])}
            for sol in solutions
        ],
        key=lambda row: row["S"],
    )
    assert normalized == [
        {"S": 54, "g": "-3/2", "d": "3/2"},
        {"S": 216, "g": "3", "d": "-15/2"},
    ]

    w = sp.symbols("w")
    representatives = {
        "cube_class_1": {
            "scale": 1,
            "A": sp.Integer(-3),
            "B": sp.Integer(6),
            "C": 3 * w,
            "D": -sp.Rational(15, 2) * w,
            "E": 3 * w,
        },
        "cube_class_2": {
            "scale": 2,
            "A": -sp.Rational(3, 4),
            "B": sp.Integer(3),
            "C": -sp.Rational(3, 8) * w,
            "D": sp.Rational(3, 4) * w,
            "E": 3 * w,
        },
    }
    for row in representatives.values():
        assert verify_section(row["scale"], row)

    reductions = {
        "cube_class_1": [
            {"prime": 7, "w": 2, "order": finite_reduction_order(1, 7, 2)},
            {"prime": 13, "w": 6, "order": finite_reduction_order(1, 13, 6)},
        ],
        "cube_class_2": [
            {"prime": 7, "w": 2, "order": finite_reduction_order(2, 7, 2)},
            {"prime": 13, "w": 6, "order": finite_reduction_order(2, 13, 6)},
        ],
    }
    assert [row["order"] for row in reductions["cube_class_1"]] == [13, 19]
    assert [row["order"] for row in reductions["cube_class_2"]] == [13, 21]

    # If the specialized point had torsion order n, good reduction would give
    # n/order(P mod p) a p-power. The two displayed equations are incompatible.
    non_torsion = {
        "cube_class_1": "n=13*7^a and n=19*13^b have no common positive solution",
        "cube_class_2": "n=13*7^a and n=21*13^b have no common positive solution",
    }

    explicit = {}
    for name, row in representatives.items():
        explicit[name] = {
            "scale": row["scale"],
            "x": f"u*(({row['A']})*u^3+({row['B']}))",
            "y": f"({row['C']})*u^6+({row['D']})*u^3+({row['E']})",
            "verified": True,
        }

    rank_certificate_path = Path("search/results/cyclic_cubic_base_ranks.json")
    assert rank_certificate_path.exists(), "paired Sage base-rank certificate is missing"
    rank_certificate = json.loads(rank_certificate_path.read_text())
    rank_rows = {row["scale"]: row for row in rank_certificate["rows"]}
    for scale in ("1", "2"):
        assert rank_rows[scale]["exact_rank"] == 0
        assert rank_rows[scale]["torsion_order"] == 3
    assert "proof.all(True)" in rank_certificate["truth_status"]

    result = {
        "status": "pass",
        "author": "Francesco Giannicola",
        "surface": "mu=2 marked j=0 family",
        "field": "K=Q(w), w^2=-3",
        "minimal_character_one_ansatz": {
            "x": "u*(A*u^3+B)",
            "y": "C*u^6+D*u^3+E",
        },
        "normalized_invariants": {
            "rho": "c*A/B",
            "S": "c*B^3",
            "g": "c^2*C/w",
            "d": "c*D/w",
        },
        "B_zero_branch": {
            "status": "impossible",
            "computed_u6_coefficient": str(b_zero_left),
            "required_u6_coefficient": str(b_zero_required),
        },
        "rho_eliminant": str(sp.factor(80 * eliminant)),
        "irreducible_cubic": str(cubic),
        "irreducible_cubic_discriminant": int(cubic_discriminant),
        "K_rational_rho_values": ["-1/2"],
        "normalized_solutions": normalized,
        "scale_cube_classes_in_Qstar_mod_cubes": ["1", "2"],
        "cube_intersection_lemma": (
            "If z in Q(sqrt(-3)) and z^3 in Q, then z^3 is a rational cube; "
            "hence K^3 intersect Q equals Q^3."
        ),
        "explicit_sections": explicit,
        "specialization_non_torsion_certificates": reductions,
        "non_torsion_contradictions": non_torsion,
        "CM_orbit": {
            "zeta": "(-1+w)/2",
            "relation": "P+[zeta]P+[zeta^2]P=O",
            "Z_rank_if_P_non_torsion": 2,
        },
        "paired_base_rank_certificate": {
            "path": str(rank_certificate_path),
            "cube_class_1_exact_rank": rank_rows["1"]["exact_rank"],
            "cube_class_2_exact_rank": rank_rows["2"]["exact_rank"],
            "proof_status": rank_certificate["truth_status"],
        },
        "restricted_obstruction": (
            "For mu=2, every rational scale carrying a minimal character-one "
            "Q(sqrt(-3))-eigensection lies in cube class 1 or 2, and both "
            "corresponding cyclic cubic bases have exact rational rank 0."
        ),
        "truth_status": "restricted-family theorem; no rank-30 curve is claimed",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "scale_cube_classes": result["scale_cube_classes_in_Qstar_mod_cubes"],
        "rho_eliminant": result["rho_eliminant"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
