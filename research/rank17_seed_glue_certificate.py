#!/usr/bin/env python3
"""Exact certificate for the rank-one seed of the X(6,79) neighbor chain.

The reconstructed positive essential lattice starts from the transparent
rank-17 lattice

    N0 = A11 (+) K6,

where K6 is the explicit determinant-79 Gram matrix below.  The norm-2 root
sublattice is A11 (+) A3 (+) A2.  This script derives, using only exact rational
arithmetic, the primitive rank-one Mordell--Weil quotient, its height 79/12,
and the component classes forced at the reducible fibres.

The final conversion from lattice data to a section intersection number uses
the standard Shioda height formula for an elliptic K3 surface.  The script does
not construct a Weierstrass model and does not choose between Kodaira I3 and IV
for the A2 fibre; that distinction is invisible to the root lattice alone.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable

Matrix = list[list[int]]
RationalMatrix = list[list[Fraction]]


def a_cartan(rank: int) -> Matrix:
    matrix = [[0 for _ in range(rank)] for _ in range(rank)]
    for index in range(rank):
        matrix[index][index] = 2
        if index + 1 < rank:
            matrix[index][index + 1] = -1
            matrix[index + 1][index] = -1
    return matrix


def block_diagonal(*blocks: Matrix) -> Matrix:
    size = sum(len(block) for block in blocks)
    result = [[0 for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                result[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return result


def transpose(matrix: list[list[object]]) -> list[list[object]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def multiply(
    left: list[list[object]], right: list[list[object]]
) -> list[list[object]]:
    if not left or not right:
        return []
    if len(left[0]) != len(right):
        raise ValueError("incompatible matrix dimensions")
    right_t = transpose(right)
    return [
        [
            sum(
                (
                    left_value * right_value
                    for left_value, right_value in zip(
                        row, column, strict=True
                    )
                ),
                0,
            )
            for column in right_t
        ]
        for row in left
    ]


def quadratic(vector: list[Fraction], gram: Matrix) -> Fraction:
    return sum(
        (
            vector[row]
            * Fraction(gram[row][column])
            * vector[column]
            for row in range(len(vector))
            for column in range(len(vector))
        ),
        Fraction(0),
    )


def bareiss_determinant(matrix: Matrix) -> int:
    size = len(matrix)
    if size == 0:
        return 1
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous_pivot = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            swap_row = next(
                (
                    row
                    for row in range(column + 1, size)
                    if work[row][column] != 0
                ),
                None,
            )
            if swap_row is None:
                return 0
            work[column], work[swap_row] = (
                work[swap_row],
                work[column],
            )
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot
                    - work[row][column] * work[column][target]
                )
                if numerator % previous_pivot:
                    raise AssertionError("Bareiss division was not exact")
                work[row][target] = numerator // previous_pivot
        previous_pivot = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def solve_exact(matrix: Matrix, right_hand_side: list[int]) -> list[Fraction]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("solve_exact requires a nonempty square matrix")
    if len(right_hand_side) != size:
        raise ValueError("right-hand side dimension mismatch")
    augmented: RationalMatrix = [
        [Fraction(value) for value in row]
        + [Fraction(right_hand_side[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot_row = next(
            (
                row
                for row in range(column, size)
                if augmented[row][column] != 0
            ),
            None,
        )
        if pivot_row is None:
            raise ArithmeticError("singular exact linear system")
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[column],
        )
        pivot = augmented[column][column]
        augmented[column] = [
            value / pivot for value in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
    return [augmented[row][-1] for row in range(size)]


def lcm_denominators(values: Iterable[Fraction]) -> int:
    result = 1
    for value in values:
        result = math.lcm(result, value.denominator)
    return result


def fraction_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_vector(values: Iterable[Fraction]) -> list[str]:
    return [fraction_string(value) for value in values]


def canonical_sha256(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()


A11 = a_cartan(11)
K6: Matrix = [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, 0],
    [0, 0, -1, 8, -1, 0],
    [0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, -1, 2],
]
N0 = block_diagonal(A11, K6)
OMITTED_INDEX = 14
ROOT_INDICES = [
    index for index in range(17) if index != OMITTED_INDEX
]
ROOT_COMPONENT_SLICES = {
    "A11": slice(0, 11),
    "A3": slice(11, 14),
    "A2": slice(14, 16),
}


def principal_submatrix(matrix: Matrix, indices: list[int]) -> Matrix:
    return [
        [matrix[row][column] for column in indices]
        for row in indices
    ]


def compute_certificate() -> dict[str, object]:
    determinant_n0 = bareiss_determinant(N0)
    determinant_k6 = bareiss_determinant(K6)
    root_gram = principal_submatrix(N0, ROOT_INDICES)
    determinant_root = bareiss_determinant(root_gram)
    expected_root_gram = block_diagonal(
        a_cartan(11), a_cartan(3), a_cartan(2)
    )

    if determinant_k6 != 79:
        raise AssertionError(
            ("unexpected determinant of K6", determinant_k6)
        )
    if determinant_n0 != 948:
        raise AssertionError(
            ("unexpected determinant of N0", determinant_n0)
        )
    if root_gram != expected_root_gram:
        raise AssertionError(
            "the coordinate root lattice is not A11+A3+A2"
        )
    if determinant_root != 144:
        raise AssertionError(
            ("unexpected root determinant", determinant_root)
        )

    omitted = [Fraction(0) for _ in range(17)]
    omitted[OMITTED_INDEX] = Fraction(1)
    pairings = [
        N0[OMITTED_INDEX][index] for index in ROOT_INDICES
    ]
    correction = solve_exact(root_gram, pairings)

    correction_full = [Fraction(0) for _ in range(17)]
    for root_position, lattice_index in enumerate(ROOT_INDICES):
        correction_full[lattice_index] = correction[root_position]
    orthogonal = [
        omitted[index] - correction_full[index]
        for index in range(17)
    ]

    orthogonality_checks = [
        sum(
            (
                orthogonal[row]
                * Fraction(N0[row][root_index])
                for row in range(17)
            ),
            Fraction(0),
        )
        for root_index in ROOT_INDICES
    ]
    if any(value != 0 for value in orthogonality_checks):
        raise AssertionError(
            (
                "projection is not orthogonal to roots",
                orthogonality_checks,
            )
        )

    component_records: dict[str, dict[str, object]] = {}
    local_sum = Fraction(0)
    expected_component_data = {
        "A11": {
            "correction": [Fraction(0)] * 11,
            "order": 1,
            "norm": Fraction(0),
            "class": "0",
            "kodaira": "I12",
            "component": "0",
        },
        "A3": {
            "correction": [
                Fraction(-1, 4),
                Fraction(-1, 2),
                Fraction(-3, 4),
            ],
            "order": 4,
            "norm": Fraction(3, 4),
            "class": "-omega_3 = omega_1 in A3*/A3",
            "kodaira": "I4",
            "component": "1 (or 3 after reversing the cycle)",
        },
        "A2": {
            "correction": [Fraction(-2, 3), Fraction(-1, 3)],
            "order": 3,
            "norm": Fraction(2, 3),
            "class": "-omega_1 = omega_2 in A2*/A2",
            "kodaira": "I3 or IV",
            "component": "2 (or 1 after reversing the diagram)",
        },
    }

    for name, component_slice in ROOT_COMPONENT_SLICES.items():
        values = correction[component_slice]
        rank = len(values)
        gram = a_cartan(rank)
        norm = quadratic(values, gram)
        order = lcm_denominators(values)
        expected = expected_component_data[name]
        if (
            values != expected["correction"]
            or order != expected["order"]
            or norm != expected["norm"]
        ):
            raise AssertionError(
                (
                    "unexpected local glue data",
                    name,
                    values,
                    order,
                    norm,
                )
            )
        local_sum += norm
        component_records[name] = {
            "rank": rank,
            "root_gram": gram,
            "pairings_of_omitted_basis_vector": pairings[
                component_slice
            ],
            "root_dual_correction_coefficients": fraction_vector(
                values
            ),
            "discriminant_class_order": order,
            "discriminant_class": expected["class"],
            "local_height_contribution": fraction_string(norm),
            "kodaira_fiber_possibilities": expected["kodaira"],
            "section_component_index": expected["component"],
        }

    complement_norm = quadratic(orthogonal, N0)
    determinant_regulator = Fraction(
        determinant_n0, determinant_root
    )
    if local_sum != Fraction(17, 12):
        raise AssertionError(
            ("unexpected local contribution sum", local_sum)
        )
    if complement_norm != Fraction(79, 12):
        raise AssertionError(
            ("unexpected orthogonal complement norm", complement_norm)
        )
    if determinant_regulator != complement_norm:
        raise AssertionError(
            (
                "determinant and projection regulators disagree",
                determinant_regulator,
                complement_norm,
            )
        )

    chi = 2
    section_zero_intersection = (
        complement_norm - Fraction(2 * chi) + local_sum
    ) / 2
    if section_zero_intersection != 2:
        raise AssertionError(
            (
                "unexpected P.O intersection",
                section_zero_intersection,
            )
        )

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "x6_79_rank17_seed_glue",
        "claim_status": (
            "CERTIFIED exact lattice computation; "
            "Weierstrass realization remains open"
        ),
        "exact_claim": (
            "For the transparent seed lattice N0=A11 direct-sum K6, "
            "the primitive root lattice is A11+A3+A2 and the rank-one "
            "orthogonal projection has square 79/12.  Its discriminant "
            "classes are zero on A11, -omega_3 on A3, and -omega_1 "
            "on A2."
        ),
        "source": {
            "neighbor_chain_seed": "A11-plus-K6-det79",
            "frozen_target_gram_sha256": (
                "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295"
            ),
            "neighbor_primes": [5, 5, 2, 2, 2, 2, 2],
            "truth_note": (
                "The seed and target hashes come from the separately "
                "replayed exact neighbor-chain artifact; this certificate "
                "recomputes only the seed glue."
            ),
        },
        "essential_lattice": {
            "description": "N0=A11 direct-sum K6",
            "rank": 17,
            "gram_matrix": N0,
            "determinant": determinant_n0,
            "K6_gram_matrix": K6,
            "K6_determinant": determinant_k6,
        },
        "root_lattice": {
            "type": "A11 + A3 + A2",
            "rank": 16,
            "coordinate_indices_zero_based": ROOT_INDICES,
            "omitted_coordinate_zero_based": OMITTED_INDEX,
            "gram_matrix": root_gram,
            "determinant": determinant_root,
            "primitive_in_N0": True,
            "primitivity_witness": (
                "The root basis consists of sixteen standard coordinate "
                "vectors; the omitted standard coordinate generates "
                "N0/R as Z."
            ),
            "predicted_section_torsion_order": 1,
        },
        "rank_one_glue": {
            "omitted_basis_vector": fraction_vector(omitted),
            "pairings_with_ordered_root_basis": pairings,
            "root_dual_correction_coefficients": fraction_vector(
                correction
            ),
            "orthogonal_projection_coefficients_in_N0_basis": (
                fraction_vector(orthogonal)
            ),
            "orthogonality_checks": fraction_vector(
                orthogonality_checks
            ),
            "components": component_records,
            "sum_local_height_contributions": fraction_string(
                local_sum
            ),
            "projection_square": fraction_string(complement_norm),
            "regulator_from_determinants": fraction_string(
                determinant_regulator
            ),
            "mordell_weil_rank": 1,
            "mordell_weil_generator_height": fraction_string(
                complement_norm
            ),
        },
        "shioda_height_consequence": {
            "formula": (
                "<P,P> = 2*chi(O_X) + 2*(P.O) "
                "- sum_v contr_v(P)"
            ),
            "elliptic_K3_chi": chi,
            "height": fraction_string(complement_norm),
            "local_contribution_sum": fraction_string(local_sum),
            "forced_intersection_P_dot_O": fraction_string(
                section_zero_intersection
            ),
            "forced_component_data": {
                "A11_fiber": "identity component",
                "A3_fiber": (
                    "nonidentity class -omega_3=omega_1; "
                    "I4 index 1 up to reversal"
                ),
                "A2_fiber": (
                    "nonidentity class -omega_1=omega_2; "
                    "index 2 up to diagram reversal"
                ),
            },
        },
        "fiber_configuration_status": {
            "lattice_forces": ["A11", "A3", "A2"],
            "kodaira_forces": ["I12", "I4"],
            "kodaira_A2_ambiguity": ["I3", "IV"],
            "semistable_candidate": "I12 + I4 + I3 + 5 I1",
            "additive_A2_candidate": "I12 + I4 + IV + 4 I1",
            "truth_note": (
                "The A2 root lattice alone cannot distinguish I3 from IV.  "
                "Any claim of a unique semistable configuration requires "
                "the missing Weierstrass model or equivalent local "
                "monodromy data."
            ),
        },
        "open_bridge": (
            "Construct an explicit elliptic K3 Weierstrass model realizing "
            "these fiber and section constraints, then transport it through "
            "the frozen 5,5,2,2,2,2,2 neighbor chain."
        ),
        "conditional_assumptions": [],
        "implementation": {
            "language": "Python standard library",
            "script": "research/rank17_seed_glue_certificate.py",
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
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
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(
            arguments.compare.read_text(encoding="utf-8")
        )
        if committed != certificate:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
        print(f"matched {arguments.compare}")

    print("seed root system: A11 + A3 + A2")
    print("rank-one height: 79/12")
    print("local contributions: 0 + 3/4 + 2/3 = 17/12")
    print("forced section intersection P.O: 2")
    print(
        "A2 Kodaira type remains I3-or-IV "
        "until a model is reconstructed"
    )
    print(
        f"certificate sha256: {certificate['certificate_sha256']}"
    )


if __name__ == "__main__":
    main()
