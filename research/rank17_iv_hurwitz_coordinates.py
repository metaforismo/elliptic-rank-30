#!/usr/bin/env python3
"""Extract intrinsic almost-Belyi coordinates from exact IV surface censuses.

For each normalized finite-field surface write

    c4=(t-1)^2 A,  c6=(t-1)^2 B,
    F=(t-1)^2 A^3,  H=F-B^2=t^4 R,
    C=2AB+3(t-1)A'B-2(t-1)AB'.

The exact surface equations imply C=c3*t^3+c4*t^4.  The fourth critical point
of phi=F/H is

    e=-c3/c4,

and its branch value is lambda=phi(e) in P^1.  These coordinates are intrinsic
to the almost-Belyi presentation and may separate components better than the
coefficient projection (r,s).

The script replays all identities, records collision divisors lambda=0,1,inf,
computes fibre multiplicities of e, lambda, and (e,lambda), and gives exact
finite-field interpolation ranks for the affine (e,lambda) points.  It makes
no characteristic-zero, irreducibility, section, or rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def mod(value: int, prime: int) -> int:
    return value % prime


def trim(values: list[int], prime: int) -> list[int]:
    result = [value % prime for value in values]
    while result and result[-1] == 0:
        result.pop()
    return result


def add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    return trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ], prime)


def scale(values: list[int], scalar: int, prime: int) -> list[int]:
    return trim([scalar * value for value in values], prime)


def subtract(left: list[int], right: list[int], prime: int) -> list[int]:
    return add(left, scale(right, -1, prime), prime)


def multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    if not left or not right:
        return []
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] = (result[i+j] + a*b) % prime
    return trim(result, prime)


def power(values: list[int], exponent: int, prime: int) -> list[int]:
    result = [1]
    base = list(values)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base, prime)
        base = multiply(base, base, prime)
        remaining //= 2
    return result


def derivative(values: list[int], prime: int) -> list[int]:
    return trim([
        index * values[index] for index in range(1, len(values))
    ], prime)


def evaluate(values: list[int], argument: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(values):
        result = (result * argument + coefficient) % prime
    return result


def divmod_poly(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    numerator = trim(dividend, prime)
    denominator = trim(divisor, prime)
    if not denominator:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(numerator) < len(denominator):
        return [], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while numerator and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index+shift] = (
                numerator[index+shift] - coefficient*value
            ) % prime
        numerator = trim(numerator, prime)
    return trim(quotient, prime), numerator


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


def matrix_rank_mod(rows: list[list[int]], prime: int) -> int:
    matrix = [[value % prime for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    pivot_row = 0
    for column in range(column_count):
        selected = next((
            row for row in range(pivot_row, row_count)
            if matrix[row][column]
        ), None)
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
                (value - factor*pivot) % prime
                for value, pivot in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def projective_lambda(numerator: int, denominator: int, prime: int):
    numerator %= prime
    denominator %= prime
    if denominator == 0:
        if numerator == 0:
            raise AssertionError("undefined branch value 0/0")
        return [1, 0], None
    value = numerator * pow(denominator, -1, prime) % prime
    return [value, 1], value


def analyse_record(record: dict[str, object], prime: int) -> dict[str, object]:
    c4 = [int(value) % prime for value in record["c4_coefficients_ascending"]]
    c6 = [int(value) % prime for value in record["c6_coefficients_ascending"]]
    factor = [1, -2, 1]
    A, remainder4 = divmod_poly(c4, factor, prime)
    B, remainder6 = divmod_poly(c6, factor, prime)
    if remainder4 or remainder6 or len(A) != 7 or len(B) != 11:
        raise AssertionError(("invalid IV factorization", record["record_sha256"]))

    term1 = scale(multiply(A, B, prime), 2, prime)
    term2 = scale(
        multiply(
            [-1, 1],
            multiply(derivative(A, prime), B, prime),
            prime,
        ),
        3,
        prime,
    )
    term3 = scale(
        multiply(
            [-1, 1],
            multiply(A, derivative(B, prime), prime),
            prime,
        ),
        -2,
        prime,
    )
    C = add(add(term1, term2, prime), term3, prime)
    support = [index for index, value in enumerate(C) if value % prime]
    if support != [3, 4]:
        raise AssertionError((
            "Wronskian bracket support is not {3,4}",
            prime,
            record["record_sha256"],
            support,
        ))
    c3 = C[3] % prime
    c4_coefficient = C[4] % prime
    if not c3 or not c4_coefficient:
        raise AssertionError("exact endpoint coefficient vanished")
    extra_point = -c3 * pow(c4_coefficient, -1, prime) % prime
    if evaluate(C, extra_point, prime) != 0:
        raise AssertionError("extra critical point does not kill C")

    A3 = power(A, 3, prime)
    F = multiply(factor, A3, prime)
    G = power(B, 2, prime)
    H = subtract(F, G, prime)
    if any(H[index] if index < len(H) else 0 for index in range(4)):
        raise AssertionError("H is not divisible by t^4")
    residual = trim(H[4:], prime)
    if len(residual) != 5:
        raise AssertionError(("residual is not quartic", residual))
    if c4_coefficient != (-12 * residual[-1]) % prime:
        raise AssertionError("Wronskian leading coefficient relation failed")

    f_value = evaluate(F, extra_point, prime)
    g_value = evaluate(G, extra_point, prime)
    h_value = evaluate(H, extra_point, prime)
    if (f_value - g_value - h_value) % prime:
        raise AssertionError("F-G-H evaluation identity failed")
    lambda_projective, lambda_affine = projective_lambda(
        f_value, h_value, prime
    )

    a_value = evaluate(A, extra_point, prime)
    b_value = evaluate(B, extra_point, prime)
    r_value = evaluate(residual, extra_point, prime)
    zero_collision = extra_point == 1 or a_value == 0
    one_collision = b_value == 0
    pole_collision = r_value == 0
    if sum((zero_collision, one_collision, pole_collision)) > 1:
        raise AssertionError("extra critical point lies on multiple branch divisors")
    if lambda_affine == 0 and not zero_collision:
        raise AssertionError("lambda=0 without a zero-fibre collision")
    if lambda_affine == 1 and not one_collision:
        raise AssertionError("lambda=1 without a one-fibre collision")
    if lambda_projective == [1, 0] and not pole_collision:
        raise AssertionError("lambda=infinity without a residual pole collision")
    if zero_collision and lambda_affine != 0:
        raise AssertionError("zero collision has nonzero branch value")
    if one_collision and lambda_affine != 1:
        raise AssertionError("one collision has branch value different from one")
    if pole_collision and lambda_projective != [1, 0]:
        raise AssertionError("pole collision has finite branch value")

    collision_type = (
        "lambda=0" if zero_collision else
        "lambda=1" if one_collision else
        "lambda=infinity" if pole_collision else
        "generic"
    )
    parameters = {
        key: int(value)
        for key, value in record["parameters"].items()
    }
    result = {
        "source_record_sha256": record["record_sha256"],
        "parameters": parameters,
        "wronskian_C_coefficients_ascending": C,
        "extra_critical_point_e": extra_point,
        "extra_branch_value_projective": lambda_projective,
        "extra_branch_value_affine": lambda_affine,
        "collision_type": collision_type,
        "A_at_e": a_value,
        "B_at_e": b_value,
        "R_at_e": r_value,
        "residual_quartic_coefficients_ascending": residual,
    }
    result["record_sha256"] = canonical_hash(result)
    return result


def fibre_summary(values: list[object]) -> dict[str, object]:
    rendered = [json.dumps(value, sort_keys=True) for value in values]
    counts = Counter(rendered)
    multiplicities = sorted(counts.values(), reverse=True)
    return {
        "source_count": len(values),
        "distinct_value_count": len(counts),
        "maximum_fibre_size": max(multiplicities, default=0),
        "fibre_size_histogram": {
            str(size): multiplicities.count(size)
            for size in sorted(set(multiplicities))
        },
    }


def affine_rank_profile(records: list[dict[str, object]], prime: int):
    points = [
        (
            int(record["extra_critical_point_e"]),
            int(record["extra_branch_value_affine"]),
        )
        for record in records
        if record["extra_branch_value_affine"] is not None
    ]
    profile = []
    for degree in range(1, 7):
        basis = monomials(degree)
        matrix = [
            [
                pow(x, i, prime) * pow(y, j, prime) % prime
                for i, j in basis
            ]
            for x, y in points
        ]
        rank = matrix_rank_mod(matrix, prime)
        profile.append({
            "total_degree": degree,
            "affine_point_count": len(points),
            "monomial_count": len(basis),
            "matrix_rank": rank,
            "nullity": len(basis) - rank,
            "sample_sufficient_for_detection": len(points) >= len(basis),
            "detected_common_relation": (
                len(points) >= len(basis) and rank < len(basis)
            ),
        })
    return {
        "finite_lambda_point_count": len(points),
        "infinite_lambda_point_count": len(records) - len(points),
        "distinct_affine_pair_count": len(set(points)),
        "rank_profile": profile,
    }


def analyse_census(path: Path) -> dict[str, object]:
    source = json.loads(path.read_text(encoding="utf-8"))
    prime = int(source["prime"])
    records = [analyse_record(record, prime) for record in source["records"]]
    e_values = [record["extra_critical_point_e"] for record in records]
    lambda_values = [
        record["extra_branch_value_projective"] for record in records
    ]
    pairs = [
        [
            record["extra_critical_point_e"],
            record["extra_branch_value_projective"],
        ]
        for record in records
    ]
    collision_counts = Counter(record["collision_type"] for record in records)
    result = {
        "prime": prime,
        "source_record_sha256": source["record_sha256"],
        "surface_count": len(records),
        "records": records,
        "extra_point_map": fibre_summary(e_values),
        "branch_value_map": fibre_summary(lambda_values),
        "hurwitz_pair_map": fibre_summary(pairs),
        "collision_counts": {
            key: collision_counts.get(key, 0)
            for key in ("generic", "lambda=0", "lambda=1", "lambda=infinity")
        },
        "affine_e_lambda_interpolation": affine_rank_profile(records, prime),
    }
    result["record_sha256"] = canonical_hash(result)
    return result


def build(paths: list[Path]) -> dict[str, object]:
    analyses = [analyse_census(path) for path in paths]
    analyses.sort(key=lambda item: item["prime"])
    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_hurwitz_coordinates",
        "truth_status": (
            "EXACT finite-field extraction of the extra critical point and "
            "branch value from certified full IV surface loci; fibre counts "
            "and interpolation ranks do not prove characteristic-zero "
            "irreducibility, rationality, a section, or rank 30"
        ),
        "primes": [item["prime"] for item in analyses],
        "analyses": analyses,
        "interpretation_boundary": [
            "Injectivity on finite-field rational points does not prove a birational map over Q.",
            "A common interpolation relation on finitely many primes is not ideal membership.",
            "Collision values 0, 1, and infinity are valid Hurwitz boundary divisors and need not be surface degenerations; residual-pole collision is separately detected.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
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
    print(json.dumps({
        "primes": payload["primes"],
        "maps": {
            str(item["prime"]): {
                "surface_count": item["surface_count"],
                "e": item["extra_point_map"],
                "lambda": item["branch_value_map"],
                "pair": item["hurwitz_pair_map"],
                "collisions": item["collision_counts"],
            }
            for item in payload["analyses"]
        },
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
