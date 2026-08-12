#!/usr/bin/env python3
"""Evaluate one exact formal-arc candidate in the F_7 rank-17 incidence germ.

The script imports the certified incidence equations from
``rank17_section_local_branch_f7.py``.  Fourteen pivot variables are solved
coefficient-by-coefficient using the invertible Jacobian minor at the certified
F_7 point.  The output records every nonzero coefficient in the five residual
equations.

This is exact arithmetic in characteristic 7.  It is not a p-adic or
characteristic-zero lifting theorem.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("rank17_section_local_branch_f7.py")


def load_module():
    spec = importlib.util.spec_from_file_location("rank17_local", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(
    a_coefficients: dict[int, int],
    w_coefficients: dict[int, int],
    order: int,
) -> tuple[dict[str, list[list[int]]], list[list[list[int]]]]:
    module = load_module()
    jacobian = module.jacobian_at_seed()
    pivot_matrix = [
        [jacobian[row][column] for column in module.PIVOT_COLUMNS]
        for row in module.PIVOT_ROWS
    ]
    variables = [module.Series(value, order, 1) for value in module.SEED]

    # a = y1-1, u=tau, v=-tau+w; hence w=u+v.
    for degree, coefficient in a_coefficients.items():
        if coefficient % 7:
            variables[14].coefficients[(degree,)] = coefficient % 7
    variables[15].coefficients[(1,)] = 1
    variables[16].coefficients[(1,)] = 6
    for degree, coefficient in w_coefficients.items():
        if coefficient % 7:
            variables[16].coefficients[(degree,)] = coefficient % 7

    pivot_coefficients: list[list[list[int]]] = []
    for degree in range(1, order + 1):
        equations = module.equations(variables)
        correction = module.solve_square(
            pivot_matrix,
            [
                -equations[row].coefficient((degree,))
                for row in module.PIVOT_ROWS
            ],
        )
        pivot_coefficients.append(
            [[column, coefficient] for column, coefficient in zip(
                module.PIVOT_COLUMNS, correction, strict=True
            ) if coefficient]
        )
        for column, coefficient in zip(
            module.PIVOT_COLUMNS, correction, strict=True
        ):
            if not coefficient:
                continue
            value = (
                variables[column].coefficient((degree,)) + coefficient
            ) % 7
            if value:
                variables[column].coefficients[(degree,)] = value
            else:
                variables[column].coefficients.pop((degree,), None)

    equations = module.equations(variables)
    residuals = {
        str(row): [
            [degree, equations[row].coefficient((degree,))]
            for degree in range(order + 1)
            if equations[row].coefficient((degree,))
        ]
        for row in module.RESIDUAL_ROWS
    }
    return residuals, pivot_coefficients


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--a", type=int, required=True)
    parser.add_argument("--w", type=int, required=True)
    parser.add_argument("--order", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    state = json.loads(arguments.state.read_text(encoding="utf-8"))
    a_coefficients = {
        int(key): int(value) % 7 for key, value in state["a"].items()
    }
    w_coefficients = {
        int(key): int(value) % 7 for key, value in state["w"].items()
    }
    a_coefficients[arguments.degree] = arguments.a % 7
    w_coefficients[arguments.degree] = arguments.w % 7

    residuals, pivots = evaluate(
        a_coefficients, w_coefficients, arguments.order
    )
    payload = {
        "schema_version": 1,
        "truth_status": (
            "VERIFIED COMPUTATION over F_7; no p-adic or "
            "characteristic-zero claim"
        ),
        "state": str(arguments.state),
        "degree": arguments.degree,
        "candidate": {
            "a": arguments.a % 7,
            "w": arguments.w % 7,
        },
        "order": arguments.order,
        "residuals": residuals,
        "first_nonzero": {
            row: values[0] if values else None
            for row, values in residuals.items()
        },
        "pivot_corrections_by_degree": pivots,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "candidate": payload["candidate"],
        "order": payload["order"],
        "first_nonzero": payload["first_nonzero"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
