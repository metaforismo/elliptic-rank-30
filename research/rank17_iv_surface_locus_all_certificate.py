#!/usr/bin/env python3
"""Certify the full normalized additive-IV surface locus over a finite field.

The production tuples come from the independent C enumerator with e0=1 and no
split-square filtering.  Every accepted tuple is rebuilt with the Python v2
invariant implementation.  A small prime is exhausted independently in Python
and must agree exactly with the C tuple set and all four split categories.

This is finite-field evidence only.  It does not prove a characteristic-zero
surface, section, or rank-30 curve.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER_PATH = HERE / "rank17_iv_surface_locus_finite_field_v2.py"


def load_helper():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_surface_locus_finite_field_v2", HELPER_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(HELPER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def surface_record_all(helper, prime: int, parameters: tuple[int, ...]):
    p0, p1, p2, p3, r, s = parameters
    if p0 % prime == 0:
        return None
    c4 = helper.c4_from_parameters(p0, p1, p2, p3, r, s, prime)
    if helper.evaluate(c4, 1, prime) or helper.derivative_at_one(c4, prime):
        raise AssertionError("the IV c4 parametrization is inconsistent")
    c6 = helper.c6_from_c4_and_constant(c4, pow(p0, 3, prime), prime)
    target_jets = [
        3 * p0 * p0 * p1,
        3 * (p0 * p0 * p2 + p0 * p1 * p1),
        3 * p0 * p0 * p3 + 6 * p0 * p1 * p2 + p1 * p1 * p1,
    ]
    if any(
        (c6[index] - target_jets[index - 1]) % prime
        for index in range(1, 4)
    ):
        return None
    if helper.evaluate(c6, 1, prime) or helper.derivative_at_one(c6, prime):
        return None
    b_at_one = helper.quadratic_coefficient_at_one(c6, prime)
    if b_at_one == 0:
        return None
    delta = helper.discriminant_numerator(c4, c6, prime)
    if any(helper.coefficient(delta, index) for index in range(4)):
        raise AssertionError("the I4 jet lost its forced vanishing")
    if helper.coefficient(delta, 4) == 0:
        return None
    if any(helper.coefficient(delta, index) for index in range(13, 25)):
        raise AssertionError("the infinity square-root recurrence lost degree")
    if helper.coefficient(delta, 12) == 0:
        return None
    i4_roots = helper.square_roots(-3 * p0, prime)
    iv_roots = helper.square_roots(-2 * b_at_one, prime)
    record = {
        "parameters": {
            "p0": p0,
            "p1": p1,
            "p2": p2,
            "p3": p3,
            "r": r,
            "s": s,
        },
        "c4_coefficients_ascending": c4,
        "c6_coefficients_ascending": c6,
        "discriminant_coefficients_ascending": delta + [0] * (25 - len(delta)),
        "exact_fibre_checks": {
            "delta_t4": helper.coefficient(delta, 4),
            "delta_t12": helper.coefficient(delta, 12),
            "c6_quadratic_coefficient_at_one": b_at_one,
        },
        "split_tangent_checks": {
            "i4_square_target": (-3 * p0) % prime,
            "i4_square_roots": i4_roots,
            "iv_square_target": (-2 * b_at_one) % prime,
            "iv_square_roots": iv_roots,
            "i4_split": bool(i4_roots),
            "iv_split": bool(iv_roots),
        },
    }
    record["record_sha256"] = canonical_hash(record)
    return record


def tuple_with_flags(record: dict[str, object]) -> tuple[int, ...]:
    parameters = record["parameters"]
    split = record["split_tangent_checks"]
    return tuple(
        int(parameters[name])
        for name in ("p0", "p1", "p2", "p3", "r", "s")
    ) + (int(split["i4_split"]), int(split["iv_split"]))


def enumerate_python(prime: int):
    helper = load_helper()
    records = []
    visited = 0
    for p0 in range(1, prime):
        for p1, p2, p3, r, s in itertools.product(range(prime), repeat=5):
            visited += 1
            record = surface_record_all(
                helper, prime, (p0, p1, p2, p3, r, s)
            )
            if record is not None:
                records.append(record)
    records.sort(key=tuple_with_flags)
    return visited, records


def category_counts(records):
    counts = {"neither": 0, "iv_only": 0, "i4_only": 0, "both": 0}
    for record in records:
        split = record["split_tangent_checks"]
        key = {
            (False, False): "neither",
            (False, True): "iv_only",
            (True, False): "i4_only",
            (True, True): "both",
        }[(bool(split["i4_split"]), bool(split["iv_split"]))]
        counts[key] += 1
    return counts


def parse_c_tuples(data):
    return {
        tuple(int(value) for value in values)
        for values in data["tuples"]
    }


def build(production_path: Path, small_path: Path):
    helper = load_helper()
    production = json.loads(production_path.read_text(encoding="utf-8"))
    small = json.loads(small_path.read_text(encoding="utf-8"))
    prime = int(production["prime"])
    small_prime = int(small["prime"])
    if production["normalization"] != "e0=1 jet-sign representative":
        raise AssertionError("unexpected production normalization")
    if small["normalization"] != production["normalization"]:
        raise AssertionError("small-prime normalization mismatch")

    production_records = []
    for values in production["tuples"]:
        tuple_values = tuple(int(value) for value in values)
        if len(tuple_values) != 8:
            raise AssertionError(("unexpected C tuple", tuple_values))
        parameters = tuple_values[:6]
        record = surface_record_all(helper, prime, parameters)
        if record is None:
            raise AssertionError(("C accepted a rejected tuple", tuple_values))
        if tuple_with_flags(record) != tuple_values:
            raise AssertionError(("C/Python flag mismatch", tuple_values))
        production_records.append(record)
    production_records.sort(key=tuple_with_flags)
    if len(production_records) != int(production["candidate_count"]):
        raise AssertionError("production candidate count mismatch")
    if category_counts(production_records) != production["split_category_counts"]:
        raise AssertionError("production category count mismatch")

    python_visited, python_small_records = enumerate_python(small_prime)
    python_small_tuples = {tuple_with_flags(record) for record in python_small_records}
    if python_visited != int(small["visited_parameter_tuples"]):
        raise AssertionError("small-prime visited count mismatch")
    if python_small_tuples != parse_c_tuples(small):
        raise AssertionError("independent small-prime tuple sets differ")
    if category_counts(python_small_records) != small["split_category_counts"]:
        raise AssertionError("independent small-prime categories differ")

    projections = {
        "r-s": sorted({
            (record["parameters"]["r"], record["parameters"]["s"])
            for record in production_records
        }),
        "p0-r": sorted({
            (record["parameters"]["p0"], record["parameters"]["r"])
            for record in production_records
        }),
        "p0-s": sorted({
            (record["parameters"]["p0"], record["parameters"]["s"])
            for record in production_records
        }),
    }
    payload = {
        "schema_version": 1,
        "certificate_id": f"rank17_iv_full_surface_locus_f{prime}",
        "truth_status": (
            f"EXHAUSTIVE normalized I12+I4+IV surface locus over F_{prime} "
            "before split-square filtering; no characteristic-zero, section, "
            "or rank-30 conclusion"
        ),
        "prime": prime,
        "normalization": production["normalization"],
        "visited_parameter_tuples": production["visited_parameter_tuples"],
        "candidate_count": len(production_records),
        "split_category_counts": category_counts(production_records),
        "records": production_records,
        "projection_point_counts": {
            name: len(points) for name, points in projections.items()
        },
        "projection_points": {
            name: [list(point) for point in points]
            for name, points in projections.items()
        },
        "independent_small_prime_cross_check": {
            "prime": small_prime,
            "visited_parameter_tuples": python_visited,
            "candidate_count": len(python_small_records),
            "tuple_sets_identical": True,
            "split_category_counts": category_counts(python_small_records),
        },
        "limitations": [
            "The residual quartic squarefreeness/coprimality open conditions are not imposed because this certificate targets the Zariski surface locus.",
            "The finite-field locus may contain points that do not lift to characteristic zero.",
            "No section equations are imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", type=Path, required=True)
    parser.add_argument("--small-prime", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.production, arguments.small_prime)
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
        "prime": payload["prime"],
        "candidate_count": payload["candidate_count"],
        "split_category_counts": payload["split_category_counts"],
        "projection_point_counts": payload["projection_point_counts"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
