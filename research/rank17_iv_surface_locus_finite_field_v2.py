#!/usr/bin/env python3
"""Independent exact verifier for the normalized I12+I4+IV surface chart.

Version 2 fixes the c4 argument of the formal c6 reconstruction and is the
only Python implementation accepted by the fast certificate workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable


def trim(values: Iterable[int]) -> list[int]:
    result = list(values)
    while result and result[-1] == 0:
        result.pop()
    return result


def add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    return trim([
        ((left[index] if index < len(left) else 0)
         + (right[index] if index < len(right) else 0)) % prime
        for index in range(size)
    ])


def subtract(left: list[int], right: list[int], prime: int) -> list[int]:
    return add(left, [(-value) % prime for value in right], prime)


def multiply(
    left: list[int],
    right: list[int],
    prime: int,
    limit: int | None = None,
) -> list[int]:
    if not left or not right:
        return []
    output_length = len(left) + len(right) - 1
    if limit is not None:
        output_length = min(output_length, limit)
    result = [0] * output_length
    for left_index, left_value in enumerate(left):
        if not left_value:
            continue
        for right_index, right_value in enumerate(right):
            target = left_index + right_index
            if target >= output_length:
                break
            if right_value:
                result[target] = (
                    result[target] + left_value * right_value
                ) % prime
    return trim(result)


def power(
    values: list[int],
    exponent: int,
    prime: int,
    limit: int | None = None,
) -> list[int]:
    result = [1]
    base = list(values)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base, prime, limit)
        remaining //= 2
        if remaining:
            base = multiply(base, base, prime, limit)
    return result


def evaluate(values: list[int], argument: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(values):
        result = (result * argument + coefficient) % prime
    return result


def derivative_at_one(values: list[int], prime: int) -> int:
    return sum(
        index * coefficient
        for index, coefficient in enumerate(values)
    ) % prime


def quadratic_coefficient_at_one(values: list[int], prime: int) -> int:
    return sum(
        (index * (index - 1) // 2) * coefficient
        for index, coefficient in enumerate(values)
    ) % prime


def square_roots(value: int, prime: int) -> list[int]:
    target = value % prime
    return [entry for entry in range(prime) if entry * entry % prime == target]


def c4_from_parameters(
    p0: int,
    p1: int,
    p2: int,
    p3: int,
    r: int,
    s: int,
    prime: int,
) -> list[int]:
    a0 = p0 * p0 % prime
    a1 = 2 * p0 * p1 % prime
    a2 = (2 * p0 * p2 + p1 * p1) % prime
    a3 = (2 * p0 * p3 + 2 * p1 * p2) % prime
    l0 = (a0 + a1 + a2 + a3 + s + 1) % prime
    l1 = (a1 + 2 * a2 + 3 * a3 + 7 * s + 8) % prime
    a4 = (l1 - 5 * l0 + r) % prime
    a5 = (4 * l0 - l1 - 2 * r) % prime
    return [a0, a1, a2, a3, a4, a5, r % prime, s % prime, 1]


def c6_from_c4_and_constant(
    c4: list[int],
    c6_constant: int,
    prime: int,
) -> list[int]:
    inverse_two = pow(2, -1, prime)
    reversed_c4 = list(reversed(c4))
    cubed = power(reversed_c4, 3, prime, limit=12)
    cubed += [0] * (12 - len(cubed))

    square_root = [0] * 12
    square_root[0] = 1
    for order in range(1, 12):
        correction = sum(
            square_root[index] * square_root[order - index]
            for index in range(1, order)
        ) % prime
        square_root[order] = (
            (cubed[order] - correction) * inverse_two
        ) % prime

    c6 = [0] * 13
    c6[0] = c6_constant % prime
    c6[12] = 1
    for order in range(1, 12):
        c6[12 - order] = square_root[order]
    return c6


def discriminant_numerator(
    c4: list[int], c6: list[int], prime: int
) -> list[int]:
    return subtract(
        power(c4, 3, prime),
        power(c6, 2, prime),
        prime,
    )


def coefficient(values: list[int], index: int) -> int:
    return values[index] if index < len(values) else 0


def surface_record(
    *,
    prime: int,
    e0: int,
    parameters: tuple[int, int, int, int, int, int],
) -> dict[str, object] | None:
    p0, p1, p2, p3, r, s = parameters
    if p0 % prime == 0:
        return None

    c4 = c4_from_parameters(p0, p1, p2, p3, r, s, prime)
    if evaluate(c4, 1, prime) or derivative_at_one(c4, prime):
        raise AssertionError("the IV c4 parametrization is inconsistent")

    c6 = c6_from_c4_and_constant(
        c4,
        e0 * pow(p0, 3, prime),
        prime,
    )
    target_jets = [
        3 * e0 * p0 * p0 * p1,
        3 * e0 * (p0 * p0 * p2 + p0 * p1 * p1),
        e0 * (
            3 * p0 * p0 * p3
            + 6 * p0 * p1 * p2
            + p1 * p1 * p1
        ),
    ]
    if any(
        (c6[index] - target_jets[index - 1]) % prime
        for index in range(1, 4)
    ):
        return None
    if evaluate(c6, 1, prime) or derivative_at_one(c6, prime):
        return None

    b_at_one = quadratic_coefficient_at_one(c6, prime)
    if not b_at_one:
        return None

    delta = discriminant_numerator(c4, c6, prime)
    if any(coefficient(delta, index) for index in range(4)):
        raise AssertionError("the I4 jet did not force order at least four")
    if coefficient(delta, 4) == 0:
        return None
    if any(coefficient(delta, index) for index in range(13, 25)):
        raise AssertionError("the I12 square-root recurrence lost a top coefficient")
    if coefficient(delta, 12) == 0:
        return None

    split_i4_value = (-3 * e0 * p0) % prime
    split_iv_value = (-2 * b_at_one) % prime
    i4_roots = square_roots(split_i4_value, prime)
    iv_roots = square_roots(split_iv_value, prime)
    if not i4_roots or not iv_roots:
        return None

    record = {
        "parameters": {
            "e0": e0 % prime,
            "p0": p0,
            "p1": p1,
            "p2": p2,
            "p3": p3,
            "r": r,
            "s": s,
        },
        "c4_coefficients_ascending": c4,
        "c6_coefficients_ascending": c6,
        "surface_equations": {
            "i4_cube_jet": [
                (c6[index] - target_jets[index - 1]) % prime
                for index in range(1, 4)
            ],
            "c4_at_one": evaluate(c4, 1, prime),
            "c4_derivative_at_one": derivative_at_one(c4, prime),
            "c6_at_one": evaluate(c6, 1, prime),
            "c6_derivative_at_one": derivative_at_one(c6, prime),
        },
        "exact_fibre_checks": {
            "delta_t4": coefficient(delta, 4),
            "delta_t12": coefficient(delta, 12),
            "c6_quadratic_coefficient_at_one": b_at_one,
        },
        "split_tangent_checks": {
            "i4_square_target": split_i4_value,
            "i4_square_roots": i4_roots,
            "iv_square_target": split_iv_value,
            "iv_square_roots": iv_roots,
        },
        "discriminant_coefficients_ascending": delta + [0] * (25 - len(delta)),
    }
    raw = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    record["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return record


def tuple_key(record: dict[str, object]) -> tuple[int, ...]:
    params = record["parameters"]
    return tuple(
        int(params[name])
        for name in ("e0", "p0", "p1", "p2", "p3", "r", "s")
    )


def enumerate_surfaces(prime: int) -> dict[str, object]:
    if prime <= 3:
        raise ValueError("prime must be greater than 3")
    candidates: list[dict[str, object]] = []
    visited = 0
    for e0 in (1, prime - 1):
        for p0 in range(1, prime):
            for p1, p2, p3, r, s in itertools.product(
                range(prime), repeat=5
            ):
                visited += 1
                record = surface_record(
                    prime=prime,
                    e0=e0,
                    parameters=(p0, p1, p2, p3, r, s),
                )
                if record is not None:
                    candidates.append(record)

    candidates.sort(key=tuple_key)
    payload: dict[str, object] = {
        "schema_version": 2,
        "certificate_id": f"rank17_iv_surface_locus_f{prime}_python_v2",
        "truth_status": (
            f"EXHAUSTIVE normalized split I12+I4+IV surface enumeration over F_{prime}; "
            "no characteristic-zero or rank-30 conclusion"
        ),
        "prime": prime,
        "visited_parameter_tuples": visited,
        "candidate_count": len(candidates),
        "candidates": candidates,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=5)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    payload = enumerate_surfaces(arguments.prime)
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({
        "prime": payload["prime"],
        "visited_parameter_tuples": payload["visited_parameter_tuples"],
        "candidate_count": payload["candidate_count"],
        "certificate_sha256": payload["certificate_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
