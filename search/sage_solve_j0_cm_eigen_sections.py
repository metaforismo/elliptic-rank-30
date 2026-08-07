#!/usr/bin/env sage-python
"""Solve the minimal character-one section scheme over Q(sqrt(-3)).

The first target is mu=2, where q(mu)=3 and the unavoidable field condition
E^2=-q(mu)^3 is exactly Eisenstein. Exact number-field solutions only are
retained; rational descent is a separate subsequent step.
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

from sage.all import PolynomialRing, QQ, QuadraticField


def q_string(value):
    return str(value)


def solve(mu):
    started = time.time()
    mu = QQ(mu)
    q0 = mu**2 - mu + 1
    r0 = mu * (mu - 1)
    if r0 == 0:
        raise ValueError("mu must not be 0 or 1")
    a = q0
    b = q0**3 / r0

    K = QuadraticField(-3, "w")
    w = K.gen()
    square_root_minus_a = (-K(a)).sqrt(all=True)
    if not square_root_minus_a:
        return {
            "status": "constant_field_too_small",
            "mu": str(mu),
            "a": str(a),
            "field": "Q(sqrt(-3))",
            "truth_status": "exact field obstruction; no rank claim",
        }

    ring = PolynomialRing(K, names=("D", "C", "B", "A"), order="lex")
    D, C, B, A = ring.gens()
    branches = []
    all_sections = []
    for root_minus_a in square_root_minus_a:
        E_value = K(a) * root_minus_a
        equations = [
            C**2 - A**3,
            2 * C * D - 3 * A**2 * B + K(a**3),
            D**2 + 2 * C * E_value - 3 * A * B**2 - K(b**2) + K(3 * a**3),
            2 * D * E_value - B**3 + K(3 * a**3),
        ]
        ideal = ring.ideal(equations)
        branch = {
            "E": str(E_value),
            "ideal_dimension": int(ideal.dimension()),
        }
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
        "a": str(a),
        "b": str(b),
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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = solve(QQ(args.mu))
    except Exception as exc:
        result = {
            "status": "error",
            "mu": args.mu,
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
