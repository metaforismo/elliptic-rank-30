#!/usr/bin/env python3
"""Enumerate all maximal totally isotropic trace codes in E8/3E8.

The output is a canonical list of 2240 row-reduced bases.  Vectors are encoded
as base-3 integers.  This is exact finite geometry and is not a rank-30 claim.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import time

# E8 Cartan/Gram matrix, reduced modulo three.
G = (
    (2, -1, 0, 0, 0, 0, 0, 0),
    (-1, 2, -1, 0, 0, 0, 0, 0),
    (0, -1, 2, -1, 0, 0, 0, 0),
    (0, 0, -1, 2, -1, 0, 0, 0),
    (0, 0, 0, -1, 2, -1, 0, -1),
    (0, 0, 0, 0, -1, 2, -1, 0),
    (0, 0, 0, 0, 0, -1, 2, 0),
    (0, 0, 0, 0, -1, 0, 0, 2),
)
G = tuple(tuple(entry % 3 for entry in row) for row in G)
INV = (0, 1, 2)


def bilinear(left, right):
    return sum(
        left[i] * G[i][j] * right[j]
        for i in range(8)
        for j in range(8)
    ) % 3


def quadratic(vector):
    # q(v)=(v,v)/2; 2^{-1}=2 in F_3.
    return 2 * bilinear(vector, vector) % 3


def normalize(vector):
    for entry in vector:
        if entry:
            inverse = INV[entry]
            return tuple(inverse * value % 3 for value in vector)
    raise ValueError("zero vector has no projective normalization")


def encode(vector):
    value = 0
    for entry in vector:
        value = 3 * value + entry
    return value


def rref(rows):
    matrix = [list(row) for row in rows if any(row)]
    if not matrix:
        return ()
    pivot_row = 0
    for column in range(8):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = INV[matrix[pivot_row][column]]
        matrix[pivot_row] = [inverse * entry % 3 for entry in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row != pivot_row and matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    (matrix[row][j] - factor * matrix[pivot_row][j]) % 3
                    for j in range(8)
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    result = [tuple(row) for row in matrix if any(row)]
    result.sort(key=lambda vector: next(i for i, entry in enumerate(vector) if entry))
    return tuple(result)


def iter_bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask -= low


def enumerate_codes():
    points = []
    for vector in product(range(3), repeat=8):
        if not any(vector) or quadratic(vector):
            continue
        if vector == normalize(vector):
            points.append(vector)
    points.sort(key=encode)
    assert len(points) == 1120

    index = {vector: i for i, vector in enumerate(points)}
    all_points = (1 << len(points)) - 1
    orthogonal_masks = []
    for left in points:
        mask = 0
        for j, right in enumerate(points):
            if bilinear(left, right) == 0:
                mask |= 1 << j
        orthogonal_masks.append(mask)

    spaces = {(point,) for point in points}
    counts = {1: len(spaces)}
    expected = {2: 36400, 3: 44800, 4: 2240}
    for target_dimension in (2, 3, 4):
        next_spaces = set()
        for basis in spaces:
            candidate_mask = all_points
            for row in basis:
                candidate_mask &= orthogonal_masks[index[row]]
            for point_index in iter_bits(candidate_mask):
                extended = rref(basis + (points[point_index],))
                if len(extended) == target_dimension:
                    next_spaces.add(extended)
        spaces = next_spaces
        counts[target_dimension] = len(spaces)
        assert len(spaces) == expected[target_dimension]

    maximal_spaces = sorted(
        spaces,
        key=lambda basis: tuple(encode(vector) for vector in basis),
    )
    encoded = [[encode(vector) for vector in basis] for basis in maximal_spaces]
    payload = json.dumps(encoded, separators=(",", ":")).encode()
    return counts, encoded, hashlib.sha256(payload).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    counts, encoded, digest = enumerate_codes()
    result = {
        "status": "pass",
        "field": 3,
        "dimension": 8,
        "witt_index": 4,
        "projective_isotropic_points": 1120,
        "subspace_counts_by_vector_dimension": {
            str(dimension): count for dimension, count in counts.items()
        },
        "maximal_space_count": 2240,
        "canonical_rref_bases_base3_encoded": encoded,
        "basis_payload_sha256": digest,
        "elapsed_seconds": time.time() - started,
        "truth_status": "exact finite geometry; no rank-30 curve claimed",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "maximal_space_count": result["maximal_space_count"],
                "basis_payload_sha256": result["basis_payload_sha256"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
