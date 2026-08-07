#!/usr/bin/env sage-python
"""Solve low-height C3-eigen polynomial sections after the cubic base change.

Start from the marked j=0 E8 rational surface

    y^2 = x^3 + b^2 r(t)^2 - a^3 q(t)^3,
    r=t(t-1), q=r+1.

After the cyclic cubic base change u^3=r, the curve descends to Q(u):

    y^2 = x^3 + b^2 u^6 - a^3 (u^3+1)^3.

For the automorphism u -> zeta_3 u, search every polynomial section disjoint
from the zero section in the three x-character classes. Exact rational
solutions only are retained.
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

from sage.all import PolynomialRing, QQ


def rational_string(value):
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def solve_character(mu, character):
    started = time.time()
    mu = QQ(mu)
    q0 = mu**2 - mu + 1
    r0 = mu * (mu - 1)
    if r0 == 0:
        raise ValueError("mu must not be 0 or 1")
    a = q0
    b = q0**3 / r0

    if character in (0, 1):
        names = ("E", "D", "C", "B", "A")
    elif character == 2:
        names = ("E", "D", "C", "A")
    else:
        raise ValueError("character must be 0, 1, or 2")

    ring = PolynomialRing(QQ, names=names, order="lex")
    variables = ring.gens_dict()
    series = PolynomialRing(ring, "u")
    u = series.gen()

    A = variables["A"]
    if character == 0:
        B = variables["B"]
        x = A * u**3 + B
    elif character == 1:
        B = variables["B"]
        x = u * (A * u**3 + B)
    else:
        x = A * u**2

    C = variables["C"]
    D = variables["D"]
    E = variables["E"]
    y = C * u**6 + D * u**3 + E
    coefficient = b**2 * u**6 - a**3 * (u**3 + 1) ** 3
    difference = y**2 - x**3 - coefficient
    equations = [ring(difference[index]) for index in range(difference.degree() + 1)]
    equations = [equation for equation in equations if equation]
    ideal = ring.ideal(equations)

    result = {
        "status": "started",
        "mu": rational_string(mu),
        "character": character,
        "a": rational_string(a),
        "b": rational_string(b),
        "curve_a6": str(coefficient),
        "x_ansatz": str(x),
        "y_ansatz": str(y),
        "equation_count": len(equations),
        "ideal_dimension": int(ideal.dimension()),
        "truth_status": "exact rational polynomial sections; independence requires a separate certificate",
    }

    if ideal.dimension() != 0:
        result["status"] = "non_zero_dimensional"
        result["elapsed_seconds"] = time.time() - started
        return result

    groebner = ideal.groebner_basis()
    result["groebner_basis_size"] = len(groebner)
    result["scheme_degree"] = int(ideal.vector_space_dimension())
    result["univariate_factors"] = []
    for polynomial in groebner:
        used = [variable for variable in ring.gens() if polynomial.degree(variable) > 0]
        if len(used) == 1:
            result["univariate_factors"].append(
                {
                    "variable": str(used[0]),
                    "degree": int(polynomial.degree()),
                    "factors": [(str(factor), int(exponent)) for factor, exponent in polynomial.factor()],
                }
            )

    try:
        variety = ideal.variety(QQ)
    except Exception as exc:
        variety = []
        result["variety_error"] = repr(exc)

    solutions = []
    section_pairs = {}
    for solution in variety:
        if any(equation.subs(solution) != 0 for equation in equations):
            raise AssertionError("Sage returned a non-solution")
        coefficients = {str(variable): rational_string(solution[variable]) for variable in ring.gens()}
        x_value = x.subs(solution)
        y_value = y.subs(solution)
        if series(y_value**2 - x_value**3 - coefficient) != 0:
            raise AssertionError("section substitution failed")
        solutions.append(
            {
                "coefficients": coefficients,
                "x": str(x_value),
                "y": str(y_value),
            }
        )
        key = str(x_value)
        section_pairs.setdefault(key, []).append(str(y_value))

    result["rational_solution_count"] = len(solutions)
    result["rational_section_pair_count"] = len(section_pairs)
    result["solutions"] = solutions
    result["section_pairs"] = [
        {"x": x_value, "y_representatives": y_values}
        for x_value, y_values in sorted(section_pairs.items())
    ]
    result["status"] = "pass"
    result["elapsed_seconds"] = time.time() - started
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", required=True)
    parser.add_argument("--character", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    try:
        result = solve_character(QQ(args.mu), args.character)
    except Exception as exc:
        result = {
            "status": "error",
            "mu": args.mu,
            "character": args.character,
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
                "character": result.get("character"),
                "scheme_degree": result.get("scheme_degree"),
                "rational_section_pairs": result.get("rational_section_pair_count"),
                "error": result.get("error"),
                "elapsed_seconds": result.get("elapsed_seconds"),
            },
            sort_keys=True,
        )
    )
    return 0 if result.get("status") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
