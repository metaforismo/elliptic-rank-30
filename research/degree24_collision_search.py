#!/usr/bin/env python3
"""Deterministic exact search for collisions in the degree-(2,4) family.

A nontrivial equality Phi(t1)=Phi(t2), where

    Phi(t)=t^12/186624-t^3/6,

produces two decompositions of the same cubic coefficient p.  The search also
checks the equivalent point on y^2=x^3-1296 and both cube-coset conditions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List


def parameter_values(numerator_bound: int, denominator_bound: int) -> Iterable[Fraction]:
    if numerator_bound < 1 or denominator_bound < 1:
        raise ValueError("bounds must be positive")
    for denominator in range(1, denominator_bound + 1):
        for numerator in range(-numerator_bound, numerator_bound + 1):
            if numerator == 0 or math.gcd(abs(numerator), denominator) != 1:
                continue
            yield Fraction(numerator, denominator)


def phi(parameter: Fraction | int) -> Fraction:
    parameter = Fraction(parameter)
    if not parameter:
        raise ValueError("parameter must be nonzero")
    return parameter**12 / 186624 - parameter**3 / 6


def rational_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def exact_cube_root(value: Fraction) -> Fraction | None:
    if not value:
        return Fraction(0)

    def integer_cube_root(number: int) -> int | None:
        sign = -1 if number < 0 else 1
        target = abs(number)
        lo, hi = 0, 1
        while hi**3 < target:
            hi *= 2
        while lo <= hi:
            mid = (lo + hi) // 2
            cube = mid**3
            if cube == target:
                return sign * mid
            if cube < target:
                lo = mid + 1
            else:
                hi = mid - 1
        return None

    numerator_root = integer_cube_root(value.numerator)
    denominator_root = integer_cube_root(value.denominator)
    if numerator_root is None or denominator_root is None:
        return None
    return Fraction(numerator_root, denominator_root)


def collision_witness(first: Fraction, second: Fraction) -> dict:
    if first == second:
        raise ValueError("collision must use distinct parameters")
    common_value = phi(first)
    if phi(second) != common_value:
        raise AssertionError("not a collision")

    q1 = first**3 / 36
    q2 = second**3 / 36
    if q1 == q2:
        raise AssertionError("distinct rational parameters gave equal q")
    relation = 3 * (q1 + q2) * (q1**2 + q2**2)
    if relation != 2:
        raise AssertionError("quartic collision relation failed")

    sum_q = q1 + q2
    if not sum_q:
        raise AssertionError("collision relation has zero sum")
    x = 12 / sum_q
    y = 36 * (q1 - q2) / sum_q
    if y**2 != x**3 - 1296:
        raise AssertionError("Mordell-curve witness failed")

    cube1 = exact_cube_root(36 * q1)
    cube2 = exact_cube_root(36 * q2)
    if cube1 != first or cube2 != second:
        raise AssertionError("cube-coset reconstruction failed")

    return {
        "mordell_point": {"x": rational_string(x), "y": rational_string(y)},
        "p": rational_string(common_value),
        "q1": rational_string(q1),
        "q2": rational_string(q2),
        "t1": rational_string(first),
        "t2": rational_string(second),
    }


def run_search(numerator_bound: int, denominator_bound: int) -> dict:
    fibres: Dict[Fraction, List[Fraction]] = defaultdict(list)
    parameters = list(parameter_values(numerator_bound, denominator_bound))
    for parameter in parameters:
        fibres[phi(parameter)].append(parameter)

    collisions: list[dict] = []
    for value in sorted(fibres, key=lambda item: (item.numerator / item.denominator, item.numerator, item.denominator)):
        fibre = sorted(fibres[value])
        if len(fibre) < 2:
            continue
        for left_index, first in enumerate(fibre):
            for second in fibre[left_index + 1 :]:
                collisions.append(collision_witness(first, second))

    result = {
        "collision_count": len(collisions),
        "collisions": collisions,
        "distinct_image_values": len(fibres),
        "domain": {
            "denominator": f"1 <= d <= {denominator_bound}",
            "form": "t=n/d in lowest terms, t!=0",
            "numerator": f"|n| <= {numerator_bound}",
        },
        "map": "Phi(t)=t^12/186624-t^3/6",
        "parameter_count": len(parameters),
        "proof_status": "exact bounded computation",
        "schema_version": 1,
        "scope_limitation": "No conclusion is made outside the stated finite rational-height box.",
    }
    canonical_without_hash = json.dumps(result, sort_keys=True, separators=(",", ":")).encode("utf-8")
    result["record_sha256"] = hashlib.sha256(canonical_without_hash).hexdigest()
    return result


def canonical_payload(result: dict) -> str:
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numerator-bound", type=int, required=True)
    parser.add_argument("--denominator-bound", type=int, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--expect-collisions", type=int)
    args = parser.parse_args()

    result = run_search(args.numerator_bound, args.denominator_bound)
    if args.expect_collisions is not None and result["collision_count"] != args.expect_collisions:
        raise SystemExit(
            f"expected {args.expect_collisions} collisions, found {result['collision_count']}"
        )
    payload = canonical_payload(result)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(
        "VERIFIED collision search",
        f"parameters={result['parameter_count']}",
        f"collisions={result['collision_count']}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
