#!/usr/bin/env sage-python
"""Solve the minimal j=0 character-one scheme over Q(sqrt(-3)) after

    u^3 = c*t*(t-1).

The positive-rank base c=36 is the primary target. Exact number-field
solutions only are retained. Rational descent through the j=0 three-isogeny
is performed by a separate verifier.
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

from sage.all import PolynomialRing, QQ, QuadraticField


def solve(mu, scale):
    started = time.time()
    mu = QQ(mu)
    scale = QQ(scale)
    if scale == 0:
        raise ValueError("scale must be nonzero")
    q0 = mu**2 - mu + 1
    r0 = mu * (mu - 1)
    if r0 == 0:
        raise ValueError("mu must not be 0 or 1")
    a = q0
    b = q0**3 / r0

    K = QuadraticField(-3, "w")
    roots = (-K(a)).sqrt(all=True)
    if not roots:
        return {
            "status": "constant_field_too_small",
            "mu": str(mu),
            "scale": str(scale),
            "a": str(a),
            "field": "Q(sqrt(-3))",
            "truth_status": "exact field obstruction; no rank claim",
        }

    ring = PolynomialRing(K, names=("D", "C", "B", "A"), order="lex")
    D, C, B, A = ring.gens()
    branches = []
    all_sections = []
    for root_minus_a in roots:
        E_value = K(a) * root_minus_a
        equations = [
            C**2 - A**3,
            2 * C * D - 3 * A**2 * B + K(a**3 / scale**3),
            D**2 + 2 * C * E_value - 3 * A * B**2 - K(b**2 / scale**2) + K(3 * a**3 / scale**2),
            2 * D * E_value - B**3 + K(3 * a**3 / scale),
        ]
        ideal = ring.ideal(equations)
        branch = {"E": str(E_value), "ideal_dimension": int(ideal.dimension())}
        if ideal.dimension() != 0:
            branch["status"] = "non_zero_dimensional"
            branches.append(branch)
            continue

        groebner = ideal.groebner_basis()
        branch["groebner_basis_size"] = len(groebner)
        branch["scheme_degree"] = int(ideal.vector_space_dimension())
        branch["univariate_factors"] = []
        for polynomial in groebner:
            used = [variable for variable in ring.gens() if polynomial.degree(variable) > 0]
            if len(used) == 1:
                branch["univariate_factors"].append(
                    {
                        "variable": str(used[0]),
                        "degree": int(polynomial.degree()),
                        "factors": [(str(factor), int(exponent)) for factor, exponent in polynomial.factor()],
                    }
                )

        try:
            variety = ideal.variety(K)
        except Exception as exc:
            variety = []
            branch["variety_error"] = repr(exc)

        sections = []
        for solution in variety:
            if any(equation.subs(solution) != 0 for equation in equations):
                raise AssertionError("number-field variety returned a non-solution")
            record = {
                "A": str(solution[A]),
                "B": str(solution[B]),
                "C": str(solution[C]),
                "D": str(solution[D]),
                "E": str(E_value),
                "x": f"u*(({solution[A]})*u^3+({solution[B]}))",
                "y": f"({solution[C]})*u^6+({solution[D]})*u^3+({E_value})",
            }
            sections.append(record)
            all_sections.append(record)
        branch["solution_count"] = len(sections)
        branch["sections"] = sections
        branch["status"] = "pass"
        branches.append(branch)

    return {
        "status": "pass",
        "mu": str(mu),
        "scale": str(scale),
        "a": str(a),
        "b": str(b),
        "base_curve": f"u^3={scale}*t*(t-1)",
        "field": "Q(sqrt(-3))",
        "branches": branches,
        "total_eigen_sections": len(all_sections),
        "sections": all_sections,
        "elapsed_seconds": time.time() - started,
        "truth_status": "exact Eisenstein eigen-sections; rational descent and independence remain open",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", default="2")
    parser.add_argument("--scale", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = solve(QQ(args.mu), QQ(args.scale))
    except Exception as exc:
        result = {
            "status": "error",
            "mu": args.mu,
            "scale": args.scale,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "truth_status": "no mathematical claim",
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": result.get("status"),
                "mu": result.get("mu"),
                "scale": result.get("scale"),
                "total_eigen_sections": result.get("total_eigen_sections"),
                "error": result.get("error"),
                "elapsed_seconds": result.get("elapsed_seconds"),
            },
            sort_keys=True,
        )
    )
    return 0 if result.get("status") in {"pass", "constant_field_too_small"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
