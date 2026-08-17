#!/usr/bin/env python3
"""Exact formal elimination at the certified local F_7 incidence point.

The 19-equation surface-plus-section system has a 14 by 14 Jacobian minor
which is invertible at the certified seed.  The formal implicit-function
theorem therefore eliminates those fourteen regular coordinates uniquely in
terms of

    a = y1 - 1,
    u = y2 - 4,
    v = y3 - 2.

This module performs that elimination degree by degree in the truncated power
series ring F_7[[a,u,v]] / (a,u,v)^(D+1).  It makes no p-adic or
characteristic-zero existence claim.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Iterator

P = 7
HERE = Path(__file__).resolve().parent
BRANCH_SCRIPT = HERE / "rank17_section_local_branch_f7.py"


def load_branch_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_section_local_branch_f7", BRANCH_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(BRANCH_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def monomials_exact(degree: int, variable_count: int = 3) -> Iterator[tuple[int, ...]]:
    if variable_count == 1:
        yield (degree,)
        return
    for first in range(degree + 1):
        for rest in monomials_exact(degree - first, variable_count - 1):
            yield (first,) + rest


def serialize_series(series, prime: int = P) -> list[dict[str, object]]:
    return [
        {"monomial": list(monomial), "coefficient": coefficient % prime}
        for monomial, coefficient in sorted(
            series.coefficients.items(), key=lambda item: (sum(item[0]), item[0])
        )
        if coefficient % prime
    ]


def compute_elimination(maximum_total_degree: int = 13) -> dict[str, object]:
    if maximum_total_degree < 1:
        raise ValueError("maximum_total_degree must be positive")

    module = load_branch_module()
    variable_count = len(module.SEED)

    dual_variables = []
    for index, value in enumerate(module.SEED):
        gradient = [0] * variable_count
        gradient[index] = 1
        dual_variables.append(module.Dual(value, gradient))
    dual_equations = module.equations(dual_variables)
    jacobian = [
        [entry % P for entry in equation.gradient]
        for equation in dual_equations
    ]
    pivot_jacobian = [
        [jacobian[row][column] for column in module.PIVOT_COLUMNS]
        for row in module.PIVOT_ROWS
    ]

    Series = module.Series
    variables = [
        Series(value, maximum_total_degree, 3) for value in module.SEED
    ]
    a = Series({(1, 0, 0): 1}, maximum_total_degree, 3)
    u = Series({(0, 1, 0): 1}, maximum_total_degree, 3)
    v = Series({(0, 0, 1): 1}, maximum_total_degree, 3)
    variables[14] = Series(module.SEED[14], maximum_total_degree, 3) + a
    variables[15] = Series(module.SEED[15], maximum_total_degree, 3) + u
    variables[16] = Series(module.SEED[16], maximum_total_degree, 3) + v

    for degree in range(1, maximum_total_degree + 1):
        equations = module.equations(variables)
        for monomial in monomials_exact(degree):
            right_hand_side = [
                -equations[row].coefficient(monomial) % P
                for row in module.PIVOT_ROWS
            ]
            coefficients = module.solve_square(
                pivot_jacobian, right_hand_side
            )
            for column, coefficient in zip(
                module.PIVOT_COLUMNS, coefficients, strict=True
            ):
                if coefficient % P:
                    variables[column].coefficients[monomial] = coefficient % P

    equations = module.equations(variables)
    pivot_remainders = []
    for row in module.PIVOT_ROWS:
        for monomial, coefficient in equations[row].coefficients.items():
            if coefficient % P:
                pivot_remainders.append(
                    {
                        "row": row,
                        "monomial": list(monomial),
                        "coefficient": coefficient % P,
                    }
                )
    if pivot_remainders:
        raise AssertionError(
            f"formal pivot equations did not vanish: {pivot_remainders[:5]}"
        )

    residual_series = []
    for row in module.RESIDUAL_ROWS:
        terms = serialize_series(equations[row])
        residual_series.append(
            {
                "row": row,
                "terms": terms,
                "minimum_degree": min(
                    (sum(term["monomial"]) for term in terms),
                    default=None,
                ),
                "term_count": len(terms),
            }
        )

    pivot_series = []
    for column in module.PIVOT_COLUMNS:
        displacement = variables[column] - Series(
            module.SEED[column], maximum_total_degree, 3
        )
        pivot_series.append(
            {
                "variable": module.VARIABLES[column],
                "terms": serialize_series(displacement),
            }
        )

    payload: dict[str, object] = {
        "schema_version": 1,
        "truth_status": (
            "CERTIFIED truncated formal elimination over F_7; "
            "no p-adic or characteristic-zero section claim"
        ),
        "prime": P,
        "maximum_total_degree": maximum_total_degree,
        "seed": module.SEED,
        "free_variables": ["a=y1-1", "u=y2-4", "v=y3-2"],
        "pivot_rows": module.PIVOT_ROWS,
        "pivot_columns": module.PIVOT_COLUMNS,
        "residual_rows": module.RESIDUAL_ROWS,
        "pivot_jacobian": pivot_jacobian,
        "pivot_series": pivot_series,
        "residual_series": residual_series,
        "limitations": [
            "Only coefficients through the declared total degree are computed.",
            "The calculation is entirely in characteristic 7.",
        ],
    }
    raw = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=13)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    payload = compute_elimination(arguments.degree)
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "status": "completed",
                "degree": arguments.degree,
                "minimum_degrees": [
                    record["minimum_degree"]
                    for record in payload["residual_series"]
                ],
                "record_sha256": payload["record_sha256"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
