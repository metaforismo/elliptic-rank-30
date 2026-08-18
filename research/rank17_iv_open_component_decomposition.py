#!/usr/bin/env python3
"""Exact determinant-open component projection for the additive-IV surface system.

The script derives the three residual coefficient polynomials at t^8,t^7,t^6
on D != 0, computes their two elimination resultants, and proves that the only
common one-dimensional projection factors are the known rational quintic F and
the determinant-zero projection G5. It also reconstructs two exact rational
zero-dimensional residual points and proves that both have identically zero
discriminant.

This does not yet classify every zero-dimensional cofactor point and therefore
does not close the whole IV branch or prove rank 30.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import subprocess
import sys
import tempfile
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


def polynomial_metadata(poly) -> dict[str, object]:
    return {
        "total_degree": int(poly.total_degree()),
        "degrees": {str(gen): int(poly.degree(gen)) for gen in poly.gens},
        "term_count": len(poly.terms()),
        "sha256": hashlib.sha256(str(poly.as_expr()).encode()).hexdigest(),
    }


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
    if remainder.as_expr() or qbar.degree() != 13:
        raise AssertionError("the exact P divisor failed")

    m2_formula = 8*(-k + 2*r3 + 4)
    m1_formula = -(
        k**2 + 14*k*r3 + 4*k + 24*r2
        + 13*r3**2 - 56*r3 - 80
    ) / 3
    q = sp.Poly(
        sp.expand(qbar.as_expr().subs({m2: m2_formula, m1: m1_formula})),
        t,
    )
    equations = [sp.expand(q.coeff_monomial(t**d)) for d in (11, 10, 9)]
    matrix, rhs = sp.linear_eq_to_matrix(equations, [m0, r1, r0])
    determinant = sp.factor(matrix.det())
    D = sp.expand(
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    if sp.cancel(determinant / D).free_symbols:
        raise AssertionError("unexpected lower determinant")
    solution = matrix.inv()*rhs
    lower = {
        m0: sp.cancel(solution[0]),
        r1: sp.cancel(solution[1]),
        r0: sp.cancel(solution[2]),
    }

    residuals = {}
    residual_denominators = {}
    for degree in (8, 7, 6):
        coefficient = sp.cancel(q.coeff_monomial(t**degree).subs(lower, simultaneous=True))
        numerator, denominator = sp.together(coefficient).as_numer_denom()
        residuals[degree] = primitive_poly(numerator, r2, r3, k)
        residual_denominators[degree] = str(sp.factor(denominator))

    P8, P7, P6 = residuals[8], residuals[7], residuals[6]
    with tempfile.TemporaryDirectory(prefix="rank17-iv-resultants-") as tmp:
        input_path = Path(tmp) / "residuals.pkl"
        output_path = Path(tmp) / "resultants.pkl"
        input_path.write_bytes(pickle.dumps(residuals, protocol=4))
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--resultant-worker",
                str(input_path),
                str(output_path),
            ],
            check=True,
        )
        resultant_87, resultant_86 = pickle.loads(output_path.read_bytes())
    resultant_87 = sp.Poly(resultant_87.as_expr(), k, r3, domain=sp.QQ)
    resultant_86 = sp.Poly(resultant_86.as_expr(), k, r3, domain=sp.QQ)
    common = sp.gcd(resultant_87, resultant_86)

    F = sp.Poly(
        81*k**5 + 1053*k**4*r3 - 1098*k**4
        + 13248*k**3*r3 - 11048*k**3
        + 9000*k**2*r3**2 + 116688*k**2*r3 + 115872*k**2
        + 12000*k*r3**2 + 87168*k*r3 + 112512*k
        + 40000*r3**3 + 384000*r3**2 + 979968*r3 + 813056,
        k, r3, domain=sp.QQ,
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
    expected_common = (F * G5**2).monic()
    if common.monic() != expected_common:
        raise AssertionError("the common resultant factor is not F*G5^2")
    common_factorization = sp.factor_list(common.as_expr())
    if sorted(
        (sp.Poly(factor, k, r3).total_degree(), exponent)
        for factor, exponent in common_factorization[1]
    ) != [(5, 1), (5, 2)]:
        raise AssertionError("unexpected common resultant factorization")

    cofactor_87 = sp.exquo(resultant_87, common)
    cofactor_86 = sp.exquo(resultant_86, common)
    if sp.gcd(cofactor_87, cofactor_86).total_degree() != 0:
        raise AssertionError("the zero-dimensional projection cofactors share a curve")

    isolated_points = []
    for coordinates in (
        {k: sp.Integer(6), r3: sp.Integer(-2), r2: sp.Integer(1)},
        {k: sp.Integer(11), r3: sp.Integer(-1), r2: sp.Integer(0)},
    ):
        lower_values = {
            m2: sp.factor(m2_formula.subs(coordinates)),
            m1: sp.factor(m1_formula.subs(coordinates)),
            m0: sp.factor(lower[m0].subs(coordinates)),
            r1: sp.factor(lower[r1].subs(coordinates)),
            r0: sp.factor(lower[r0].subs(coordinates)),
        }
        all_values = {**coordinates, **lower_values}
        nonzero_coefficients = {
            degree: str(sp.factor(q.coeff_monomial(t**degree).subs(all_values)))
            for degree in range(14)
            if sp.factor(q.coeff_monomial(t**degree).subs(all_values)) != 0
        }
        if nonzero_coefficients:
            raise AssertionError(("isolated point missed a residual", nonzero_coefficients))

        R_point = sp.expand(R.subs(all_values))
        P_point = sp.expand(P.subs(all_values))
        Q_point = sp.expand(Q.subs(all_values))
        K_point = sp.expand(K.subs(all_values))
        M_point = sp.expand(M.subs(all_values))
        N_point = sp.expand(Q_point**2 + M_point*P_point)
        T_point = sp.expand(
            N_point*K_point*(s*Q_point - 2*P_point)
            - 3*P_point*s*(
                sp.diff(N_point, t)*K_point
                - 2*N_point*sp.diff(K_point, t)
            )
        )
        A_point, rem_A = sp.div(sp.Poly(N_point, t), sp.Poly(K_point**2, t))
        B_point, rem_B = sp.div(sp.Poly(T_point, t), sp.Poly(K_point**4, t))
        if rem_A.as_expr() or rem_B.as_expr():
            raise AssertionError("isolated point failed A/B divisibility")
        A_expr = sp.factor(A_point.as_expr())
        B_expr = sp.factor(B_point.as_expr())
        H = sp.factor(s**2*A_expr**3 - B_expr**2)
        if H != 0:
            raise AssertionError("the isolated residual point is not discriminant-zero")
        isolated_points.append({
            "coordinates": {
                "k": str(coordinates[k]),
                "r3": str(coordinates[r3]),
                "r2": str(coordinates[r2]),
            },
            "lower_coordinates": {
                "m2": str(lower_values[m2]),
                "m1": str(lower_values[m1]),
                "m0": str(lower_values[m0]),
                "r1": str(lower_values[r1]),
                "r0": str(lower_values[r0]),
            },
            "D": str(sp.factor(D.subs(coordinates))),
            "A": str(A_expr),
            "B": str(B_expr),
            "discriminant_numerator": "0",
            "classification": "DEGENERATE_IDENTICALLY_ZERO_DISCRIMINANT",
        })

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_determinant_open_component_projection",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "truth_status": (
            "EXACT projection-level decomposition of three determinant-open "
            "differential residuals. It identifies every common curve factor, "
            "but does not yet classify all zero-dimensional cofactor points; no "
            "section or rank-30 conclusion."
        ),
        "determinant": str(D),
        "residuals": {
            str(degree): {
                **polynomial_metadata(poly),
                "denominator": residual_denominators[degree],
            }
            for degree, poly in residuals.items()
        },
        "resultants": {
            "Res_r2(P8,P7)": polynomial_metadata(resultant_87),
            "Res_r2(P8,P6)": polynomial_metadata(resultant_86),
            "common_factor": str(common.monic().as_expr()),
            "common_factorization": "F * G5^2",
            "F": str(F.as_expr()),
            "G5": str(G5.as_expr()),
            "cofactor_gcd": "1",
            "cofactor_degrees": {
                "P8_P7": int(cofactor_87.total_degree()),
                "P8_P6": int(cofactor_86.total_degree()),
            },
        },
        "isolated_rational_residual_points": isolated_points,
        "mathematical_consequence": (
            "The only one-dimensional projection components common to the t8, "
            "t7, and t6 residual equations are the rational quintic F and the "
            "determinant-zero branch G5. Any further determinant-open solution "
            "is zero-dimensional in projection. Two exact rational residual "
            "points in that cofactor scheme are reconstructed, and both have "
            "identically zero discriminant."
        ),
        "highest_value_next_action": (
            "Compute a rational-univariate representation of the zero-dimensional "
            "cofactor scheme after saturating by F*G5, then prove that its only "
            "Q-rational points are the two discriminant-zero records above."
        ),
        "limitations": [
            "A gcd of plane resultants classifies common curve factors, not all isolated intersections.",
            "The determinant-zero G5 branch is excluded by a separate exact rational obstruction.",
            "The rational F component is excluded from the split seed search by a separate real square-class obstruction.",
            "No height-79/12 section is constructed.",
            "The unconditional global rank lower bound remains 29.",
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def markdown(payload: dict[str, object]) -> str:
    return "\n".join([
        "# Determinant-open additive-IV projection checkpoint", "",
        "```text", "SOLVED-A: false", "SOLVED-B: false",
        "unconditional global lower bound: rank E(Q) >= 29", "```", "",
        payload["truth_status"], "", "## Consequence", "",
        payload["mathematical_consequence"], "", "## Highest-value next action", "",
        payload["highest_value_next_action"], "", "## Exact record", "",
        "```json", json.dumps(payload, indent=2, sort_keys=True), "```", "",
    ])


def main() -> None:
    if len(sys.argv) == 4 and sys.argv[1] == "--resultant-worker":
        source = pickle.loads(Path(sys.argv[2]).read_bytes())
        p8, p7, p6 = source[8], source[7], source[6]
        out = []
        elimination_variable = p8.gens[0]
        remaining_generators = p8.gens[1:]
        for other in (p7, p6):
            raw = sp.resultant(p8.as_expr(), other.as_expr(), elimination_variable)
            result = sp.Poly(raw, *remaining_generators, domain=sp.QQ)
            _, result = result.clear_denoms(convert=True)
            _, result = result.primitive()
            if result.LC() < 0:
                result = -result
            out.append(result)
        Path(sys.argv[3]).write_bytes(pickle.dumps(tuple(out), protocol=4))
        return
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
        "common_factorization": payload["resultants"]["common_factorization"],
        "cofactor_gcd": payload["resultants"]["cofactor_gcd"],
        "isolated_rational_point_count": len(payload["isolated_rational_residual_points"]),
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
