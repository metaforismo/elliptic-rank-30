#!/usr/bin/env python3
"""Exact first-order lift obstruction for the target I4/I3 component seed.

The 19 surface-plus-section coefficient equations alone are nonreduced at the
certified F_7 seed.  They do not encode, as radical equations, that the section
passes through the two multiplicative nodes and follows the required tangent
branches.  This certificate adds those six exact local equations:

    X(0) = -p0,             X(1) = -q0,
    Y(0) = 0,               Y(1) = 0,
    Y'(0) = X'(0)+p1,       Y'(1) = X'(1)+q1.

The signs in the two tangent equations are the certified component ratios
(+1,+1) of the modular seed.  For the resulting 25 equations in 17 variables,
the Jacobian over F_7 has full column rank 17.  Nevertheless the unique
first-order compatibility system for a lift modulo 7^2 is inconsistent.  An
explicit left-kernel vector ell satisfies

    ell * J = 0,      ell * rhs = 3 != 0  (mod 7).

Consequently there is no solution modulo 49 reducing to this exact seed and
satisfying these target-component equations.  Hence there is no Z_7 point in
this component/sign chart above the seed.  This does not exclude other modular
seeds, other primes, the IV branch, or a rank-30 curve.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Callable

P = 7
ROOT = Path(__file__).resolve().parents[1]
BRANCH_SCRIPT = ROOT / "research" / "rank17_section_local_branch_f7.py"
VARIABLES = [
    "p0", "p1", "p2", "p3", "q0", "q1", "q2", "s",
    "x0", "x1", "x2", "x3", "x4", "y0", "y1", "y2", "y3",
]
EQUATION_LABELS = (
    [f"surface_{index}" for index in range(6)]
    + [f"section_coefficient_{index}" for index in range(13)]
    + [
        "node_X_at_0",
        "node_X_at_1",
        "node_Y_at_0",
        "node_Y_at_1",
        "tangent_ratio_at_0_plus_1",
        "tangent_ratio_at_1_plus_1",
    ]
)


def load_branch_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_section_local_branch_f7", BRANCH_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(BRANCH_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def convolution(left: list[int], right: list[int], length: int, modulus: int) -> list[int]:
    result = [0] * length
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            if left_index + right_index < length:
                result[left_index + right_index] = (
                    result[left_index + right_index]
                    + left_value * right_value
                ) % modulus
    return result


def coefficient_equations(values: list[int], modulus: int) -> list[int]:
    """Evaluate the original 19 equations modulo an odd modulus."""

    inverse_two = pow(2, -1, modulus)
    p0, p1, p2, p3, q0, q1, q2, s, *rest = [
        value % modulus for value in values
    ]
    x_values = rest[:5] + [0] * 4
    y_values = rest[5:9] + [0] * 9

    c4 = [0] * 9
    c4[0] = p0 * p0 % modulus
    c4[1] = 2 * p0 * p1 % modulus
    c4[2] = (2 * p0 * p2 + p1 * p1) % modulus
    c4[3] = (2 * p0 * p3 + 2 * p1 * p2) % modulus
    c4[4] = (
        -3 - 15 * p0 * p0 - 20 * p0 * p1 - 12 * p0 * p2
        - 6 * p0 * p3 - 6 * p1 * p1 - 6 * p1 * p2
        + 15 * q0 * q0 - 10 * q0 * q1 + 2 * q0 * q2
        + q1 * q1 - s
    ) % modulus
    c4[5] = (
        8 + 24 * p0 * p0 + 30 * p0 * p1 + 16 * p0 * p2
        + 6 * p0 * p3 + 8 * p1 * p1 + 6 * p1 * p2
        - 24 * q0 * q0 + 18 * q0 * q1 - 4 * q0 * q2
        - 2 * q1 * q1 + 3 * s
    ) % modulus
    c4[6] = (
        -6 - 10 * p0 * p0 - 12 * p0 * p1 - 6 * p0 * p2
        - 2 * p0 * p3 - 3 * p1 * p1 - 2 * p1 * p2
        + 10 * q0 * q0 - 8 * q0 * q1 + 2 * q0 * q2
        + q1 * q1 - 3 * s
    ) % modulus
    c4[7] = s
    c4[8] = 1

    reversed_c4 = list(reversed(c4))
    c4_cubed = convolution(
        convolution(reversed_c4, reversed_c4, 12, modulus),
        reversed_c4,
        12,
        modulus,
    )
    square_root = [1]
    for order in range(1, 12):
        correction = sum(
            square_root[index] * square_root[order - index]
            for index in range(1, order)
        ) % modulus
        square_root.append(
            (c4_cubed[order] - correction) * inverse_two % modulus
        )

    c6 = [0] * 13
    c6[0] = p0**3 % modulus
    c6[12] = 1
    for order in range(1, 12):
        c6[12 - order] = square_root[order]

    result = [
        c6[1] - 3 * p0 * p0 * p1,
        c6[2] - 3 * (p0 * p0 * p2 + p0 * p1 * p1),
        c6[3] - (3 * p0 * p0 * p3 + 6 * p0 * p1 * p2 + p1**3),
        sum(c6) - q0**3,
        sum(index * c6[index] for index in range(13))
        - 3 * q0 * q0 * q1,
        sum(
            (index * (index - 1) // 2) * c6[index]
            for index in range(13)
        ) - 3 * (q0 * q0 * q2 + q0 * q1 * q1),
    ]

    y_squared = convolution(y_values, y_values, 13, modulus)
    x_squared = convolution(x_values, x_values, 13, modulus)
    x_cubed = convolution(x_squared, x_values, 13, modulus)
    c4_times_x = convolution(c4, x_values, 13, modulus)
    result.extend(
        y_squared[index]
        - x_cubed[index]
        + 3 * c4_times_x[index]
        + 2 * c6[index]
        for index in range(13)
    )
    return [value % modulus for value in result]


def component_equations(values: list[int], modulus: int) -> list[int]:
    """Radical node and signed tangent equations for D=1."""

    (
        p0, p1, _p2, _p3, q0, q1, _q2, _s,
        x0, x1, x2, x3, x4, y0, y1, y2, y3,
    ) = [value % modulus for value in values]
    return [
        x0 + p0,
        x0 + x1 + x2 + x3 + x4 + q0,
        y0,
        y0 + y1 + y2 + y3,
        y1 - x1 - p1,
        y1 + 2 * y2 + 3 * y3
        - x1 - 2 * x2 - 3 * x3 - 4 * x4 - q1,
    ]


def full_equations(values: list[int], modulus: int) -> list[int]:
    return [
        value % modulus
        for value in (
            coefficient_equations(values, modulus)
            + component_equations(values, modulus)
        )
    ]


def jacobian_mod_p(
    function: Callable[[list[int], int], list[int]],
    seed: list[int],
) -> list[list[int]]:
    """Recover the exact Jacobian mod p by first differences mod p^2."""

    modulus = P * P
    base = function(seed, modulus)
    matrix = [[0] * len(seed) for _ in base]
    for column in range(len(seed)):
        displaced = list(seed)
        displaced[column] += P
        values = function(displaced, modulus)
        for row in range(len(base)):
            difference = (values[row] - base[row]) % modulus
            if difference % P:
                raise AssertionError("first difference is not divisible by p")
            matrix[row][column] = difference // P % P
    return matrix


def rref(matrix: list[list[int]]) -> tuple[list[list[int]], list[int]]:
    work = [[value % P for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    pivots: list[int] = []
    row = 0
    for column in range(column_count):
        selected = next(
            (index for index in range(row, row_count) if work[index][column]),
            None,
        )
        if selected is None:
            continue
        work[row], work[selected] = work[selected], work[row]
        inverse = pow(work[row][column], -1, P)
        work[row] = [value * inverse % P for value in work[row]]
        for other in range(row_count):
            if other != row and work[other][column]:
                factor = work[other][column]
                work[other] = [
                    (value - factor * pivot_value) % P
                    for value, pivot_value in zip(
                        work[other], work[row], strict=True
                    )
                ]
        pivots.append(column)
        row += 1
        if row == row_count:
            break
    return work, pivots


def nullspace(matrix: list[list[int]]) -> list[list[int]]:
    reduced, pivots = rref(matrix)
    column_count = len(matrix[0]) if matrix else 0
    free_columns = [
        column for column in range(column_count) if column not in pivots
    ]
    basis = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free] % P
        basis.append(vector)
    return basis


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def dot(left: list[int], right: list[int]) -> int:
    return sum(a * b for a, b in zip(left, right, strict=True)) % P


def vector_times_matrix(
    vector: list[int], matrix: list[list[int]]
) -> list[int]:
    return [
        sum(vector[row] * matrix[row][column] for row in range(len(matrix)))
        % P
        for column in range(len(matrix[0]))
    ]


def canonical_hash(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def compute_certificate() -> dict[str, object]:
    module = load_branch_module()
    seed = list(module.SEED)
    if seed != [2, 2, 1, 0, 2, 2, 2, 0, 5, 6, 2, 0, 6, 0, 1, 4, 2]:
        raise AssertionError("unexpected modular seed")
    if any(full_equations(seed, P)):
        raise AssertionError("the target-component seed is not on the F7 system")

    jacobian = jacobian_mod_p(full_equations, seed)
    _reduced, pivot_columns = rref(jacobian)
    rank = len(pivot_columns)
    if rank != 17 or pivot_columns != list(range(17)):
        raise AssertionError(("unexpected full-system Jacobian rank", rank, pivot_columns))

    values_mod_49 = full_equations(seed, P * P)
    if any(value % P for value in values_mod_49):
        raise AssertionError("seed values modulo 49 are not divisible by 7")
    rhs = [-(value // P) % P for value in values_mod_49]

    left_kernel = nullspace(transpose(jacobian))
    if len(left_kernel) != len(EQUATION_LABELS) - rank:
        raise AssertionError("unexpected left-kernel dimension")
    incompatible = [vector for vector in left_kernel if dot(vector, rhs)]
    if len(incompatible) != 1:
        raise AssertionError(("unexpected incompatible witness count", len(incompatible)))
    witness = incompatible[0]
    witness_times_jacobian = vector_times_matrix(witness, jacobian)
    witness_dot_rhs = dot(witness, rhs)
    if any(witness_times_jacobian) or witness_dot_rhs != 3:
        raise AssertionError("left-kernel obstruction identity failed")

    sparse_witness = [
        {
            "equation": EQUATION_LABELS[index],
            "coefficient_mod_7": coefficient,
        }
        for index, coefficient in enumerate(witness)
        if coefficient
    ]

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_target_component_mod49_obstruction",
        "truth_status": (
            "CERTIFIED no lift modulo 49, hence no Z_7 lift, of this exact "
            "F_7 seed in the required split-I4/split-I3 component-sign chart; "
            "no global rank-30 conclusion"
        ),
        "prime": P,
        "seed": {
            "variables": VARIABLES,
            "values_mod_7": seed,
            "surface_and_section_equations": 19,
            "radical_component_equations": 6,
            "component_tangent_ratios_mod_7": [1, 1],
        },
        "radical_component_equations": [
            "X(0)+p0=0",
            "X(1)+q0=0",
            "Y(0)=0",
            "Y(1)=0",
            "Y'(0)-X'(0)-p1=0",
            "Y'(1)-X'(1)-q1=0",
        ],
        "linearized_lift_problem": {
            "equation_labels": EQUATION_LABELS,
            "jacobian_shape": [25, 17],
            "jacobian_mod_7": jacobian,
            "jacobian_rank_mod_7": rank,
            "pivot_columns": pivot_columns,
            "negative_seed_residual_divided_by_7_mod_7": rhs,
            "left_kernel_dimension": len(left_kernel),
            "obstruction_witness": witness,
            "obstruction_witness_sparse": sparse_witness,
            "witness_times_jacobian_mod_7": witness_times_jacobian,
            "witness_dot_rhs_mod_7": witness_dot_rhs,
        },
        "mathematical_consequence": {
            "mod_49_lift_exists": False,
            "Z_7_lift_reducing_to_seed_exists": False,
            "reason": (
                "Any lift x=seed+7*d mod 49 would require J*d=rhs mod 7. "
                "Multiplication by the certified left-kernel witness gives "
                "0=3 mod 7, a contradiction."
            ),
            "scope": (
                "Only this exact modular seed and its certified (+1,+1) "
                "target-component chart are excluded."
            ),
        },
        "correction_to_prior_local_claim": {
            "withdrawn_claim": (
                "The 19 coefficient equations alone were previously described "
                "as having m-primary initial ideal containing t^12."
            ),
            "correction": (
                "That nilpotence artifact was never successfully generated. "
                "The strategic obstruction instead comes from adding the six "
                "radical node/tangent equations and checking mixed-characteristic "
                "first-order compatibility."
            ),
            "nonreduced_system_note": (
                "The coefficient identity alone records node incidence through "
                "powers and can exhibit long nilpotent/formal directions.  Such "
                "directions are not lifts in the required component chart."
            ),
        },
        "limitations": [
            "Other F_7 section-incidence seeds are not excluded.",
            "Other primes and the additive IV realization are not excluded.",
            "This certificate does not construct the rank-17 K3 over Q.",
            "The unconditional global lower bound remains rank at least 29.",
        ],
        "implementation": {
            "language": "Python standard library",
            "script": "research/rank17_component_lift_obstruction_f7.py",
        },
    }
    payload["certificate_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != certificate:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
    print(
        json.dumps(
            {
                "truth_status": certificate["truth_status"],
                "jacobian_rank_mod_7": certificate["linearized_lift_problem"][
                    "jacobian_rank_mod_7"
                ],
                "witness_dot_rhs_mod_7": certificate["linearized_lift_problem"][
                    "witness_dot_rhs_mod_7"
                ],
                "certificate_sha256": certificate["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
