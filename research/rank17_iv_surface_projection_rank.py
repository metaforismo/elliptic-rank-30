#!/usr/bin/env python3
"""Exact evaluation-matrix ranks for projections of the full IV surface locus.

For each exact finite-field census and each stored two-coordinate projection,
this script evaluates all monomials of total degree at most d and computes the
matrix rank over F_p.  Full column rank proves that no nonzero polynomial in
that monomial space vanishes on the entire recorded point set.

This does not prove irreducibility and does not exclude a lower-degree equation
for one proper component of a reducible locus.  No characteristic-zero,
section, or rank-30 conclusion is inferred.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def monomials(total_degree: int) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(total_degree + 1)
        for j in range(total_degree + 1 - i)
    ]


def rank_mod(rows: list[list[int]], prime: int) -> int:
    matrix = [[value % prime for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (
                row
                for row in range(pivot_row, row_count)
                if matrix[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = (
            matrix[selected], matrix[pivot_row]
        )
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [
            value * inverse % prime for value in matrix[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (value - factor * pivot) % prime
                for value, pivot in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def analyse_projection(
    points: list[list[int]], prime: int, maximum_degree: int
) -> dict[str, object]:
    normalized = [
        tuple(int(value) % prime for value in point)
        for point in points
    ]
    if any(len(point) != 2 for point in normalized):
        raise ValueError("projection points must have two coordinates")
    if len(set(normalized)) != len(normalized):
        raise AssertionError("projection point list contains duplicates")

    profile = []
    for degree in range(1, maximum_degree + 1):
        basis = monomials(degree)
        rows = [
            [
                pow(x, i, prime) * pow(y, j, prime) % prime
                for i, j in basis
            ]
            for x, y in normalized
        ]
        rank = rank_mod(rows, prime)
        nullity = len(basis) - rank
        profile.append({
            "total_degree": degree,
            "point_count": len(normalized),
            "monomial_count": len(basis),
            "matrix_rank": rank,
            "nullity": nullity,
            "full_column_rank": rank == len(basis),
            "sample_sufficient_for_detection": len(normalized) >= len(basis),
            "no_common_relation_of_degree_at_most_d": rank == len(basis),
        })
    return {
        "point_count": len(normalized),
        "rank_profile": profile,
    }


def analyse_census(path: Path, maximum_degree: int) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    prime = int(data["prime"])
    projections = {}
    for name in sorted(data["projection_points"]):
        points = data["projection_points"][name]
        if len(points) != int(data["projection_point_counts"][name]):
            raise AssertionError(("projection count mismatch", path, name))
        projections[name] = analyse_projection(
            points, prime, maximum_degree
        )
    result = {
        "prime": prime,
        "source_record_sha256": data["record_sha256"],
        "surface_count": int(data["candidate_count"]),
        "maximum_total_degree": maximum_degree,
        "projections": projections,
    }
    result["record_sha256"] = canonical_hash(result)
    return result


def build(paths: list[Path], maximum_degree: int) -> dict[str, object]:
    analyses = [analyse_census(path, maximum_degree) for path in paths]
    analyses.sort(key=lambda item: item["prime"])
    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_surface_projection_evaluation_ranks",
        "truth_status": (
            "EXACT finite-field evaluation-matrix ranks for full-locus "
            "coordinate projections; full rank excludes only a polynomial "
            "containing every recorded component at that prime and does not "
            "prove irreducibility, a characteristic-zero obstruction, a "
            "section, or rank 30"
        ),
        "maximum_total_degree": maximum_degree,
        "analyses": analyses,
        "logical_boundary": [
            "Full column rank means no nonzero polynomial in the chosen monomial space vanishes on the entire recorded point set.",
            "It does not exclude a lower-degree equation for one proper component of a reducible locus.",
            "A modular relation at one prime is not an equation over Q.",
            "No section equations are imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, action="append", required=True)
    parser.add_argument("--maximum-degree", type=int, default=6)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.census, arguments.maximum_degree)
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != payload:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
    print(json.dumps({
        "primes": [analysis["prime"] for analysis in payload["analyses"]],
        "profiles": {
            str(analysis["prime"]): analysis["projections"]
            for analysis in payload["analyses"]
        },
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
