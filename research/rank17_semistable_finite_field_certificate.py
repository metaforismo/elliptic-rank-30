#!/usr/bin/env python3
"""Recompute the finite-field seed search for the rank-17 K3 bridge.

The companion C11 program exhaustively enumerates a normalized Hermite--Padé
model for elliptic K3 invariants with discriminant orders 4, 3, and 12 at
0, 1, and infinity.  This wrapper compiles the enumerator, runs it over F_5
and F_7, performs independent exact checks on every reported example, and
emits a canonical machine-readable research certificate.

This is a VERIFIED COMPUTATION about the normalized fibre locus.  It does not
construct the missing characteristic-zero Weierstrass model, its rank-one
section, or the final rank-17 fibration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

DEFAULT_SOURCE = Path("research/rank17_semistable_fiber_search.c")
PRIMES = (5, 7)


def canonical_sha256(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def polynomial_trim(values: Iterable[int], prime: int) -> list[int]:
    result = [value % prime for value in values]
    while result and result[-1] == 0:
        result.pop()
    return result


def polynomial_multiply(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index]
                + left_value * right_value
            ) % prime
    return result


def polynomial_subtract(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (
            (left[index] if index < len(left) else 0)
            - (right[index] if index < len(right) else 0)
        ) % prime
    return polynomial_trim(result, prime)


def polynomial_divmod(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    remainder = polynomial_trim(dividend, prime)
    divisor = polynomial_trim(divisor, prime)
    if not divisor:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(remainder) < len(divisor):
        return [], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse_leading = pow(divisor[-1], -1, prime)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse_leading % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % prime
        remainder = polynomial_trim(remainder, prime)
    return polynomial_trim(quotient, prime), remainder


def polynomial_gcd(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    left = polynomial_trim(left, prime)
    right = polynomial_trim(right, prime)
    while right:
        _, remainder = polynomial_divmod(left, right, prime)
        left, right = right, remainder
    if not left:
        return []
    inverse_leading = pow(left[-1], -1, prime)
    return [(value * inverse_leading) % prime for value in left]


def taylor_jet(polynomial: list[int], order: int, prime: int) -> int:
    return sum(
        math.comb(index, order) * coefficient
        for index, coefficient in enumerate(polynomial)
        if index >= order
    ) % prime


def verify_example(experiment: dict, sign_record: dict) -> None:
    prime = experiment["prime"]
    example = sign_record["example_smooth_semistable"]
    if example is None:
        raise AssertionError("missing smooth semistable example")
    parameters = example["parameters"]
    if len(parameters) != 8:
        raise AssertionError("unexpected parameter count")
    p0, p1, p2, p3, q0, q1, q2, _s = parameters
    if p0 % prime == 0 or q0 % prime == 0:
        raise AssertionError("example violates nonvanishing constraints")

    c4 = example["c4_coefficients_ascending"]
    c6 = example["c6_coefficients_ascending"]
    if len(c4) != 9 or c4[-1] % prime != 1:
        raise AssertionError("c4 is not monic of degree 8")
    if len(c6) != 13 or c6[-1] % prime != 1:
        raise AssertionError("c6 is not monic of degree 12")

    cube_c4 = polynomial_multiply(
        polynomial_multiply(c4, c4, prime), c4, prime
    )
    square_c6 = polynomial_multiply(c6, c6, prime)
    discriminant_numerator = polynomial_subtract(
        cube_c4, square_c6, prime
    )
    if len(discriminant_numerator) != 13:
        raise AssertionError("discriminant does not have degree 12")
    if any(discriminant_numerator[index] for index in range(4)):
        raise AssertionError("discriminant order at t=0 is below 4")
    if discriminant_numerator[4] == 0:
        raise AssertionError("discriminant order at t=0 exceeds 4")
    if any(taylor_jet(discriminant_numerator, order, prime) for order in range(3)):
        raise AssertionError("discriminant order at t=1 is below 3")
    if taylor_jet(discriminant_numerator, 3, prime) == 0:
        raise AssertionError("discriminant order at t=1 exceeds 3")

    divisor = [0] * 4 + [(-1) % prime, 3 % prime, (-3) % prime, 1]
    residual, remainder = polynomial_divmod(
        discriminant_numerator, divisor, prime
    )
    if remainder:
        raise AssertionError("discriminant did not divide by t^4(t-1)^3")
    expected_residual = [
        value % prime
        for value in example[
            "residual_quintic_coefficients_ascending"
        ]
    ]
    if residual != expected_residual or len(residual) != 6:
        raise AssertionError("reported residual quintic is inconsistent")
    derivative = [
        index * residual[index] % prime
        for index in range(1, len(residual))
    ]
    if len(polynomial_gcd(residual, derivative, prime)) != 1:
        raise AssertionError("reported residual quintic is not squarefree")

    local0 = [p0, p1, p2, p3]
    local1 = [q0, q1, q2]
    local0_cube = polynomial_multiply(
        polynomial_multiply(local0, local0, prime), local0, prime
    )
    local1_cube = polynomial_multiply(
        polynomial_multiply(local1, local1, prime), local1, prime
    )
    e0 = sign_record["e0"] % prime
    e1 = sign_record["e1"] % prime
    if [value % prime for value in c6[:4]] != [
        e0 * value % prime for value in local0_cube[:4]
    ]:
        raise AssertionError("c6 does not match the signed cubic jet at 0")
    c6_at_one = [taylor_jet(c6, order, prime) for order in range(3)]
    if c6_at_one != [e1 * value % prime for value in local1_cube[:3]]:
        raise AssertionError("c6 does not match the signed cubic jet at 1")
    if example["formal_jacobian_rank"] != 6:
        raise AssertionError("example is not a rank-6 Jacobian point")


def validate_experiment(experiment: dict) -> None:
    prime = experiment["prime"]
    if prime not in PRIMES:
        raise AssertionError(f"unexpected prime: {prime}")
    expected_raw = (prime - 1) ** 2 * prime**6
    if experiment["raw_parameter_tuples"] != expected_raw:
        raise AssertionError("raw tuple count is inconsistent")
    if experiment["fiber_orders"] != {
        "t=0": 4,
        "t=1": 3,
        "t=infinity": 12,
    }:
        raise AssertionError("unexpected normalized fibre orders")
    records = experiment["sign_data"]
    if [(record["e0"], record["e1"]) for record in records] != [
        (1, 1),
        (1, -1),
        (-1, 1),
        (-1, -1),
    ]:
        raise AssertionError("sign records are not canonical")

    count_vectors = []
    for record in records:
        exact = record["exact_fixed_order_solutions"]
        rank_counts = record["formal_jacobian_rank_counts"]
        gcd_counts = record["residual_gcd_degree_counts"]
        if len(rank_counts) != 7 or sum(rank_counts) != exact:
            raise AssertionError("Jacobian-rank counts do not sum to exact count")
        if rank_counts != [0, 0, 0, 0, 0, 0, exact]:
            raise AssertionError("not every exact point has formal Jacobian rank 6")
        if len(gcd_counts) != 6 or sum(gcd_counts) != exact:
            raise AssertionError("residual gcd counts do not sum to exact count")
        if record["smooth_rank6_solutions"] != exact:
            raise AssertionError("smooth count is inconsistent")
        if record["squarefree_semistable_solutions"] != gcd_counts[0]:
            raise AssertionError("squarefree count is inconsistent")
        if (
            record["smooth_squarefree_semistable_solutions"]
            != gcd_counts[0]
        ):
            raise AssertionError("smooth squarefree count is inconsistent")
        verify_example(experiment, record)
        count_vectors.append(
            (
                record["jet_solutions"],
                exact,
                tuple(rank_counts),
                tuple(gcd_counts),
            )
        )
    if len(set(count_vectors)) != 1:
        raise AssertionError("the four sign branches have asymmetric counts")


def run_enumerator(source: Path, prime: int, executable: Path) -> dict:
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    completed = subprocess.run(
        [str(executable), str(prime)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
    )
    if completed.stderr.strip():
        raise AssertionError(
            f"enumerator wrote unexpected stderr for p={prime}: "
            f"{completed.stderr}"
        )
    result = json.loads(completed.stdout)
    validate_experiment(result)
    return result


def compute_certificate(source: Path) -> dict[str, object]:
    source = source.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    with tempfile.TemporaryDirectory(prefix="rank17-finite-field-") as directory:
        executable = Path(directory) / "rank17_semistable_fiber_search"
        subprocess.run(
            [
                "cc",
                "-O3",
                "-std=c11",
                "-Wall",
                "-Wextra",
                "-pedantic",
                str(source),
                "-o",
                str(executable),
            ],
            check=True,
        )
        experiments = [
            run_enumerator(source, prime, executable) for prime in PRIMES
        ]

    aggregate = {
        str(experiment["prime"]): {
            "exact_fixed_order_sign_parameter_tuples": sum(
                record["exact_fixed_order_solutions"]
                for record in experiment["sign_data"]
            ),
            "smooth_squarefree_semistable_sign_parameter_tuples": sum(
                record["smooth_squarefree_semistable_solutions"]
                for record in experiment["sign_data"]
            ),
        }
        for experiment in experiments
    }
    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_semistable_finite_field_locus",
        "truth_status": "VERIFIED COMPUTATION",
        "exact_claim": (
            "The normalized six-equation Hermite--Pade fibre system was "
            "exhaustively enumerated over F_5 and F_7 under p0*q0 != 0.  "
            "The committed counts include exact discriminant orders "
            "(4,3,12), formal Jacobian ranks, and squarefreeness of the "
            "residual quintic."
        ),
        "normalization": {
            "base_coordinate": "t",
            "fixed_places": ["0", "1", "infinity"],
            "target_discriminant_orders": [4, 3, 12],
            "invariants": {
                "c4": "monic polynomial of degree 8",
                "c6": "monic polynomial of degree 12",
                "discriminant_numerator": "c4^3-c6^2",
            },
            "local_square_root_jets": {
                "at_0": "l0=p0+p1*t+p2*t^2+p3*t^3",
                "at_1": "l1=q0+q1*(t-1)+q2*(t-1)^2",
                "constraints": [
                    "c4 == l0^2 mod t^4",
                    "c4 == l1^2 mod (t-1)^3",
                    "c6 == e0*l0^3 mod t^4",
                    "c6 == e1*l1^3 mod (t-1)^3",
                ],
            },
            "parameters": [
                "p0",
                "p1",
                "p2",
                "p3",
                "q0",
                "q1",
                "q2",
                "s",
            ],
            "equation_count": 6,
            "ambient_dimension": 8,
            "expected_smooth_dimension": 2,
            "sign_branches": [[1, 1], [1, -1], [-1, 1], [-1, -1]],
        },
        "experiments": experiments,
        "aggregate_counts": aggregate,
        "interpretation": {
            "jacobian_conclusion": (
                "Every exact fixed-order point found over F_5 and F_7 has "
                "formal Jacobian rank 6, so the six-equation locus is smooth "
                "of dimension 2 at each of those points."
            ),
            "semistable_conclusion": (
                "For each sign branch, 16 of 18 exact F_5 points and 41 of "
                "49 exact F_7 points have squarefree residual quintic, hence "
                "the normalized discriminant has five additional simple "
                "roots over an algebraic closure."
            ),
            "next_use": (
                "Use a smooth F_7 point as a Hensel seed, fix two transverse "
                "parameters, lift the remaining six variables p-adically, "
                "and test rational reconstruction before imposing the "
                "rank-one section glue."
            ),
        },
        "limitations": [
            (
                "Counts are sign-parameter tuples, not isomorphism classes "
                "of elliptic K3 surfaces."
            ),
            (
                "Smooth finite-field points give p-adic deformation seeds "
                "but do not by themselves produce a rational point of the "
                "characteristic-zero parameter surface."
            ),
            (
                "The experiment imposes the semistable I3 branch for the A2 "
                "root lattice; the alternative IV branch remains separate."
            ),
            (
                "No Mordell--Weil section or rank-17 Weierstrass fibration is "
                "constructed by this certificate."
            ),
        ],
        "implementation": {
            "enumerator": "research/rank17_semistable_fiber_search.c",
            "wrapper": "research/rank17_semistable_finite_field_certificate.py",
            "language": "C11 plus Python standard library",
            "compile_command": (
                "cc -O3 -std=c11 -Wall -Wextra -pedantic "
                "research/rank17_semistable_fiber_search.c"
            ),
            "enumerator_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "primes": list(PRIMES),
        },
        "conditional_assumptions": [],
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    certificate = compute_certificate(arguments.source)
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text())
        if committed != certificate:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
        print(f"matched {arguments.compare}")

    for experiment in certificate["experiments"]:
        first = experiment["sign_data"][0]
        print(
            f"F_{experiment['prime']}: exact={first['exact_fixed_order_solutions']} "
            f"per sign, smooth semistable="
            f"{first['smooth_squarefree_semistable_solutions']} per sign"
        )
    print(f"certificate sha256: {certificate['certificate_sha256']}")
    print("truth status: VERIFIED COMPUTATION; no rank-30 claim")


if __name__ == "__main__":
    main()
