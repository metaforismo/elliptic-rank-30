#!/usr/bin/env sage -python
"""Interpolate and cross-validate a plane equation for the IV surface locus.

The script uses the (r,s)-projections of exact full finite-field censuses.  A
relation of total degree d is called detected only when the number of sampled
points is at least the number of monomials and the evaluation matrix has a
nontrivial kernel.  The unique quartic through the 14 F19 points is therefore
recorded as under-sampled, while an F23 quartic is eligible for detection only
when at least 15 distinct projected points are available.

When normalized quartics at F19 and F23 exist, their coefficients are combined
by CRT, balanced in (-M/2,M/2], and tested against every projected point over
F11, F13, F17, F19, and F23.  Even complete agreement is an interpolation
hypothesis, not ideal membership over Q.  Promotion requires exact reduction
by the Wronskian elimination ideal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ
from sage.version import version as sage_version


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


def rref_nullspace_mod(
    rows: list[list[int]], prime: int
) -> tuple[int, list[list[int]]]:
    matrix = [[value % prime for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    pivot_columns = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (
                row for row in range(pivot_row, row_count)
                if matrix[row][column] % prime
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
                (value - factor * pivot_value) % prime
                for value, pivot_value in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    free_columns = [
        column for column in range(column_count)
        if column not in pivot_columns
    ]
    basis = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -matrix[row][free] % prime
        basis.append(vector)
    return len(pivot_columns), basis


def normalize_vector(
    vector: list[int], prime: int, monomial_basis: list[tuple[int, int]]
) -> tuple[list[int], tuple[int, int]]:
    preferred = (max(i for i, _j in monomial_basis), 0)
    preferred_index = (
        monomial_basis.index(preferred)
        if preferred in monomial_basis else None
    )
    if preferred_index is not None and vector[preferred_index] % prime:
        index = preferred_index
    else:
        index = max(
            position for position, value in enumerate(vector)
            if value % prime
        )
    inverse = pow(vector[index] % prime, -1, prime)
    return (
        [value * inverse % prime for value in vector],
        monomial_basis[index],
    )


def evaluate_vector(
    vector: list[int],
    monomial_basis: list[tuple[int, int]],
    point: tuple[int, int],
    prime: int,
) -> int:
    x, y = point
    return sum(
        coefficient * pow(x, i, prime) * pow(y, j, prime)
        for coefficient, (i, j) in zip(
            vector, monomial_basis, strict=True
        )
    ) % prime


def load_census(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    prime = int(data["prime"])
    points = [
        (int(point[0]), int(point[1]))
        for point in data["projection_points"]["r-s"]
    ]
    if len(points) != int(data["projection_point_counts"]["r-s"]):
        raise AssertionError(("projection count mismatch", path))
    if len(set(points)) != len(points):
        raise AssertionError(("duplicate r-s projection point", path))
    return {
        "path": str(path),
        "prime": prime,
        "record_sha256": data["record_sha256"],
        "points": points,
    }


def degree_analysis(census, maximum_degree: int = 6):
    prime = census["prime"]
    points = census["points"]
    records = []
    for degree in range(1, maximum_degree + 1):
        basis = monomials(degree)
        rows = [
            [
                pow(x, i, prime) * pow(y, j, prime) % prime
                for i, j in basis
            ]
            for x, y in points
        ]
        rank, kernel = rref_nullspace_mod(rows, prime)
        normalized = []
        for vector in kernel:
            value, pivot_monomial = normalize_vector(vector, prime, basis)
            if any(evaluate_vector(value, basis, point, prime) for point in points):
                raise AssertionError("nullspace vector failed point replay")
            normalized.append({
                "coefficients": value,
                "normalizing_monomial": list(pivot_monomial),
            })
        records.append({
            "total_degree": degree,
            "point_count": len(points),
            "monomial_count": len(basis),
            "matrix_rank": rank,
            "nullity": len(kernel),
            "sample_sufficient_for_detection": len(points) >= len(basis),
            "detected_relation": (
                len(points) >= len(basis) and len(kernel) > 0
            ),
            "monomials": [list(value) for value in basis],
            "normalized_kernel_basis": normalized,
        })
    return records


def quartic_vector(analysis):
    record = next(item for item in analysis if item["total_degree"] == 4)
    if record["nullity"] != 1:
        return None
    item = record["normalized_kernel_basis"][0]
    return {
        "coefficients": item["coefficients"],
        "monomials": record["monomials"],
        "normalizing_monomial": item["normalizing_monomial"],
        "sample_sufficient_for_detection": record[
            "sample_sufficient_for_detection"
        ],
    }


def crt_pair(left: int, p: int, right: int, q: int) -> int:
    return (
        left + p * (((right-left) * pow(p, -1, q)) % q)
    ) % (p*q)


def balanced(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if residue > modulus // 2 else residue


def vector_to_polynomial(vector, monomial_basis):
    ring = PolynomialRing(QQ, names=("r", "s"), order="lex")
    r, s = ring.gens()
    polynomial = ring(sum(
        coefficient * r**i * s**j
        for coefficient, (i, j) in zip(
            vector, monomial_basis, strict=True
        )
    ))
    content = math.gcd(*[abs(int(value)) for value in vector if value])
    if content > 1:
        polynomial = ring(polynomial / content)
        vector = [value // content for value in vector]
    if vector[-1] < 0:
        polynomial = -polynomial
        vector = [-value for value in vector]
    return ring, polynomial, vector


def build(paths: list[Path]) -> dict[str, object]:
    censuses = [load_census(path) for path in paths]
    censuses.sort(key=lambda item: item["prime"])
    by_prime = {item["prime"]: item for item in censuses}
    required = {11, 13, 17, 19, 23}
    if set(by_prime) != required:
        raise ValueError(("expected exactly F11,F13,F17,F19,F23", sorted(by_prime)))

    analyses = {
        prime: degree_analysis(census)
        for prime, census in by_prime.items()
    }
    quartics = {
        prime: quartic_vector(analyses[prime])
        for prime in sorted(analyses)
    }
    q19 = quartics[19]
    q23 = quartics[23]
    candidate = None
    if q19 is not None and q23 is not None:
        monomial_basis = [tuple(value) for value in q19["monomials"]]
        if q23["monomials"] != q19["monomials"]:
            raise AssertionError("quartic monomial bases differ")
        if q19["normalizing_monomial"] != q23["normalizing_monomial"]:
            raise AssertionError("quartic normalizations differ")
        modulus = 19 * 23
        crt_coefficients = [
            crt_pair(left, 19, right, 23)
            for left, right in zip(
                q19["coefficients"], q23["coefficients"], strict=True
            )
        ]
        integer_coefficients = [
            balanced(value, modulus) for value in crt_coefficients
        ]
        ring, polynomial, integer_coefficients = vector_to_polynomial(
            integer_coefficients, monomial_basis
        )
        reductions = {}
        all_vanish = True
        for prime, census in by_prime.items():
            values = [
                evaluate_vector(
                    [value % prime for value in integer_coefficients],
                    monomial_basis,
                    point,
                    prime,
                )
                for point in census["points"]
            ]
            vanish = all(value == 0 for value in values)
            all_vanish = all_vanish and vanish
            reductions[str(prime)] = {
                "point_count": len(values),
                "all_points_vanish": vanish,
                "nonzero_value_count": sum(value != 0 for value in values),
                "values": values,
            }
        factorization = [
            {"factor": str(factor), "exponent": int(exponent)}
            for factor, exponent in polynomial.factor()
        ]
        candidate = {
            "construction_primes": [19, 23],
            "modulus": modulus,
            "normalizing_monomial": q19["normalizing_monomial"],
            "monomials": [list(value) for value in monomial_basis],
            "balanced_integer_coefficients": integer_coefficients,
            "polynomial": str(polynomial),
            "factorization_over_Q": factorization,
            "validation_reductions": reductions,
            "all_listed_finite_field_points_vanish": all_vanish,
            "classification": (
                "MULTIPRIME_INTERPOLATION_HYPOTHESIS"
                if all_vanish else
                "CRT_INTERPOLANT_REJECTED_BY_INDEPENDENT_PRIME"
            ),
            "not_ideal_membership": True,
        }

    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_surface_plane_interpolation",
        "truth_status": (
            "Exact finite-field interpolation and cross-prime point replay; "
            "a surviving polynomial is only a characteristic-zero hypothesis "
            "until exact Wronskian-ideal membership is proved; no section or "
            "rank-30 conclusion"
        ),
        "sage_version": str(sage_version),
        "projection": ["r", "s"],
        "sources": {
            str(prime): {
                "record_sha256": census["record_sha256"],
                "point_count": len(census["points"]),
            }
            for prime, census in by_prime.items()
        },
        "degree_analyses": {
            str(prime): analyses[prime]
            for prime in sorted(analyses)
        },
        "quartic_records": {
            str(prime): quartics[prime]
            for prime in sorted(quartics)
        },
        "integer_candidate": candidate,
        "promotion_gate": (
            "Reduce the candidate polynomial to zero in the saturated "
            "Wronskian surface ideal over Q, and prove that the projection "
            "degree accounts for the intended component."
        ),
        "limitations": [
            "Point interpolation can agree at finitely many primes without giving ideal membership over Q.",
            "The F19 quartic is under-sampled because 14 points are used for 15 monomials.",
            "Bad-reduction components may be absent from individual finite-field point sets.",
            "No Mordell-Weil section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.census)
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
    candidate = payload["integer_candidate"]
    print(json.dumps({
        "candidate_classification": None if candidate is None else candidate["classification"],
        "candidate_polynomial": None if candidate is None else candidate["polynomial"],
        "all_points_vanish": None if candidate is None else candidate["all_listed_finite_field_points_vanish"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
