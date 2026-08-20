#!/usr/bin/env python3
"""Verify the frozen additive-IV component-classification certificate.

Default mode checks every compact exact identity, the determinant-zero
rational obstruction, the rational-component split obstruction, and the
Hurwitz passport.  ``--heavy`` additionally recomputes the three large
pairwise resultants of e8 with e7,e6,e5 and verifies that their common curve
factor is exactly F*S^2.

The heavy mode is intentionally separate because SymPy's exact multivariate
resultants are expensive.  The repository's Sage/Singular elimination script
is the preferred independent heavy replay.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def canonical_hash_without_record(payload: dict[str, object]) -> str:
    raw = json.dumps(
        {key: value for key, value in payload.items()
         if key != "record_sha256"},
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(raw).hexdigest()


def primitive_poly(expression, *generators) -> sp.Poly:
    poly = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    _, poly = poly.clear_denoms(convert=True)
    _, poly = poly.primitive()
    if poly.LC() < 0:
        poly = -poly
    return poly


def parse_poly(record: dict[str, object], symbols: dict[str, sp.Symbol],
               generators: tuple[sp.Symbol, ...]) -> sp.Poly:
    expression = sp.sympify(record["expression"], locals=symbols)
    poly = sp.Poly(expression, *generators, domain=sp.QQ)
    if int(record["total_degree"]) != poly.total_degree():
        raise AssertionError("stored total degree mismatch")
    if int(record["term_count"]) != len(poly.terms()):
        raise AssertionError("stored term count mismatch")
    if record["sha256"] != hashlib.sha256(
        str(poly.as_expr()).encode()
    ).hexdigest():
        raise AssertionError("stored polynomial hash mismatch")
    return poly


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--certificate",
        type=Path,
        default=Path(
            "research/data/"
            "RANK17_IV_COMPONENT_CLASSIFICATION_FROZEN.json"
        ),
    )
    parser.add_argument("--heavy", action="store_true")
    arguments = parser.parse_args()

    payload = json.loads(arguments.certificate.read_text(encoding="utf-8"))
    if canonical_hash_without_record(payload) != payload["record_sha256"]:
        raise AssertionError("global record hash mismatch")

    k, r3, r2, m0, w = sp.symbols("k r3 r2 m0 w", real=True)
    symbols = {"k": k, "r3": r3, "r2": r2, "m0": m0, "w": w}

    components = payload["components"]
    F = parse_poly(components["F"], symbols, (k, r3))
    S = parse_poly(components["S"], symbols, (k, r3))
    L = parse_poly(
        components["Res_e8_D"]["L"], symbols, (k, r3)
    )
    R1 = parse_poly(components["R1"], symbols, (r2, k, r3))
    R2 = parse_poly(components["R2"], symbols, (r2, k, r3))

    for name, polynomial in (("F", F), ("S", S)):
        factors = sp.factor_list(polynomial.as_expr())[1]
        if len(factors) != 1 or factors[0][1] != 1:
            raise AssertionError(f"{name} is reducible over Q")
        if sp.Poly(factors[0][0], k, r3).total_degree() != 5:
            raise AssertionError(f"{name} factor has wrong degree")

    resultant_e8_D = sp.Poly(
        sp.sympify(
            components["Res_e8_D"]["expression"], locals=symbols
        ),
        k, r3, domain=sp.QQ,
    )
    if primitive_poly(resultant_e8_D.as_expr(), k, r3) != primitive_poly(
        L.as_expr()*S.as_expr(), k, r3
    ):
        raise AssertionError("Res(e8,D) != L*S")

    singular = payload["singular_D_zero_chart"]
    C = parse_poly(singular["C"], symbols, (k, r3))
    G = parse_poly(singular["G"], symbols, (k, r3))
    res_SL = sp.Poly(
        sp.resultant(S.as_expr(), L.as_expr(), k),
        r3, domain=sp.QQ,
    )
    res_SC = sp.Poly(
        sp.resultant(S.as_expr(), C.as_expr(), k),
        r3, domain=sp.QQ,
    )
    if res_SL.gcd(res_SC).degree() != 0:
        raise AssertionError("special L=0 singular chart is nonempty")

    res_SG = primitive_poly(
        sp.resultant(S.as_expr(), G.as_expr(), k), r3
    )
    expected_factors = []
    for record in singular["Res_S_G"]["factors"]:
        factor = sp.Poly(
            sp.sympify(record["expression"], locals=symbols),
            r3, domain=sp.QQ,
        )
        if factor.degree() != int(record["degree"]):
            raise AssertionError("stored S-G factor degree mismatch")
        if record["sha256"] != hashlib.sha256(
            str(factor.as_expr()).encode()
        ).hexdigest():
            raise AssertionError("stored S-G factor hash mismatch")
        if sp.factor_list(factor.as_expr())[1] != [
            (factor.as_expr(), 1)
        ]:
            # SymPy can normalize the returned factor by a rational scalar.
            normalized = primitive_poly(factor.as_expr(), r3)
            listed = sp.factor_list(normalized.as_expr())[1]
            if len(listed) != 1 or listed[0][1] != 1:
                raise AssertionError("stored S-G factor is reducible")
        expected_factors.append(factor.as_expr()**int(record["exponent"]))
    expected_res_SG = primitive_poly(
        sp.prod(expected_factors), r3
    )
    if res_SG != expected_res_SG:
        raise AssertionError("Res_k(S,G) factorization mismatch")
    if any(
        int(record["degree"]) == 1
        for record in singular["Res_S_G"]["factors"]
    ):
        raise AssertionError("singular chart has a rational r3 candidate")

    rational = payload["rational_F_component"]
    k_w = sp.sympify(rational["k_w"], locals=symbols)
    r3_w = sp.sympify(rational["r3_w"], locals=symbols)
    r2_w = sp.sympify(rational["r2_w"], locals=symbols)
    if sp.cancel(F.as_expr().subs({k: k_w, r3: r3_w})) != 0:
        raise AssertionError("F parametrization failed")
    if sp.cancel(
        R1.as_expr().subs({k: k_w, r3: r3_w, r2: r2_w})
    ) != 0:
        raise AssertionError("R1 parametrization failed")
    if sp.cancel(
        R2.as_expr().subs({k: k_w, r3: r3_w, r2: r2_w})
    ) != 0:
        raise AssertionError("R2 parametrization failed")

    i4_class = sp.sympify(
        rational["I4_square_class"], locals=symbols
    )
    iv_class = sp.sympify(
        rational["IV_square_class"], locals=symbols
    )
    i4_positive = sp.solve_univariate_inequality(
        i4_class > 0, w, relational=False
    )
    iv_positive = sp.solve_univariate_inequality(
        iv_class > 0, w, relational=False
    )
    if sp.Intersection(i4_positive, iv_positive) != sp.EmptySet:
        raise AssertionError("split positivity regions overlap")

    passport = payload["hurwitz_passport"]
    partitions = [
        passport["branch_0"],
        passport["branch_1"],
        passport["branch_infinity"],
        passport["moving_branch"],
    ]
    if any(sum(partition) != 20 for partition in partitions):
        raise AssertionError("passport degree mismatch")
    ramification = sum(
        sum(index - 1 for index in partition)
        for partition in partitions
    )
    if ramification != 38 or ramification != 2*20 - 2:
        raise AssertionError("Riemann-Hurwitz mismatch")

    heavy_summary = None
    if arguments.heavy:
        residuals = {
            int(degree): parse_poly(record, symbols, (r2, r3, k))
            for degree, record in payload["residuals"].items()
        }
        e8 = residuals[8].as_expr()
        resultants = []
        for degree in (7, 6, 5):
            resultant = primitive_poly(
                sp.resultant(e8, residuals[degree].as_expr(), r2),
                k, r3,
            )
            stored = payload["projected_resultants"][
                f"Res_e8_e{degree}"
            ]
            if hashlib.sha256(
                str(resultant.as_expr()).encode()
            ).hexdigest() != stored["sha256"]:
                raise AssertionError(
                    f"heavy resultant hash mismatch at e{degree}"
                )
            resultants.append(resultant)
        common = resultants[0].gcd(resultants[1]).gcd(resultants[2])
        _, common = common.primitive()
        if common.LC() < 0:
            common = -common
        if common != F*S**2:
            raise AssertionError("heavy common factor is not F*S^2")
        heavy_summary = {
            "resultant_degrees": [
                int(resultant.total_degree())
                for resultant in resultants
            ],
            "common_factor": "F*S^2",
        }

    summary = {
        "certificate_id": payload["certificate_id"],
        "record_sha256": payload["record_sha256"],
        "F_irreducible": True,
        "S_irreducible": True,
        "singular_chart_has_Q_point": False,
        "split_positive_intersection": "EmptySet",
        "passport_ramification_total": ramification,
        "heavy_replay": heavy_summary,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
