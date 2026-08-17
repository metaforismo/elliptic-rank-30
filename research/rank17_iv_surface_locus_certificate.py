#!/usr/bin/env python3
"""Build and independently verify the finite-field IV surface certificate.

The exhaustive loop is performed by rank17_iv_surface_locus_finite_field.c.
This verifier reconstructs every accepted model using the independent Python
implementation and cross-checks the complete accepted tuple set at a smaller
prime where a full Python enumeration is inexpensive.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "rank17_iv_surface_locus_finite_field.py"


def load_implementation():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_surface_locus_finite_field", IMPLEMENTATION
    )
    if spec is None or spec.loader is None:
        raise ImportError(IMPLEMENTATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tuple_key(record: dict[str, object]) -> tuple[int, ...]:
    params = record["parameters"]
    return tuple(
        int(params[name])
        for name in ("e0", "p0", "p1", "p2", "p3", "r", "s")
    )


def expected_visited(prime: int) -> int:
    return 2 * (prime - 1) * prime**5


def load_raw(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    prime = int(data["prime"])
    tuples = [tuple(int(value) for value in item) for item in data["tuples"]]
    if int(data["visited_parameter_tuples"]) != expected_visited(prime):
        raise AssertionError("the C loop did not visit the complete chart")
    if int(data["candidate_count"]) != len(tuples):
        raise AssertionError("candidate count does not match tuple list")
    if len(tuples) != len(set(tuples)):
        raise AssertionError("duplicate C candidates")
    if tuples != sorted(tuples):
        raise AssertionError("C candidates are not in canonical order")
    return {
        "prime": prime,
        "visited_parameter_tuples": int(data["visited_parameter_tuples"]),
        "tuples": tuples,
    }


def rebuild_records(raw: dict[str, object], implementation) -> list[dict[str, object]]:
    prime = int(raw["prime"])
    records = []
    for values in raw["tuples"]:
        e0, p0, p1, p2, p3, r, s = values
        record = implementation.surface_record(
            prime=prime,
            e0=e0,
            parameters=(p0, p1, p2, p3, r, s),
        )
        if record is None:
            raise AssertionError(("C emitted a rejected tuple", values))
        if tuple_key(record) != values:
            raise AssertionError(("tuple reconstruction changed coordinates", values))
        records.append(record)
    return records


def cross_check_small_prime(
    raw_small: dict[str, object], implementation
) -> dict[str, object]:
    prime = int(raw_small["prime"])
    exhaustive = implementation.enumerate_surfaces(prime)
    python_tuples = [tuple_key(record) for record in exhaustive["candidates"]]
    c_tuples = list(raw_small["tuples"])
    if c_tuples != python_tuples:
        missing_from_c = sorted(set(python_tuples) - set(c_tuples))
        extra_in_c = sorted(set(c_tuples) - set(python_tuples))
        raise AssertionError({
            "prime": prime,
            "missing_from_c": missing_from_c[:20],
            "extra_in_c": extra_in_c[:20],
        })
    return {
        "prime": prime,
        "visited_parameter_tuples": exhaustive["visited_parameter_tuples"],
        "candidate_count": exhaustive["candidate_count"],
        "python_certificate_sha256": exhaustive["certificate_sha256"],
        "accepted_tuple_sets_identical": True,
    }


def build_certificate(raw_production: Path, raw_small: Path) -> dict[str, object]:
    implementation = load_implementation()
    production = load_raw(raw_production)
    small = load_raw(raw_small)
    if production["prime"] == small["prime"]:
        raise ValueError("production and cross-check primes must differ")

    records = rebuild_records(production, implementation)
    cross_check = cross_check_small_prime(small, implementation)
    prime = int(production["prime"])
    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": f"rank17_iv_surface_locus_f{prime}",
        "truth_status": (
            f"EXHAUSTIVE normalized split I12+I4+IV surface enumeration over F_{prime}; "
            "no characteristic-zero, rational-section, or rank-30 conclusion"
        ),
        "prime": prime,
        "visited_parameter_tuples": production["visited_parameter_tuples"],
        "candidate_count": len(records),
        "normalization": {
            "fixed_fibres": {"I4": "t=0", "IV": "t=1", "I12": "t=infinity"},
            "c4_degree": 8,
            "c6_degree": 12,
            "leading_c4": 1,
            "leading_c6": 1,
            "surface_parameters": ["e0", "p0", "p1", "p2", "p3", "r", "s"],
        },
        "independent_small_prime_cross_check": cross_check,
        "candidates": records,
        "limitations": [
            "Only the declared affine normalization and both e0 signs are covered.",
            "No section variables are imposed by this certificate.",
            "Finite-field points may fail to lift, or may lift only over number fields.",
            "The global unconditional rank lower bound remains 29.",
        ],
        "implementation": {
            "enumerator": "research/rank17_iv_surface_locus_finite_field.c",
            "independent_verifier": "research/rank17_iv_surface_locus_finite_field.py",
            "certificate_builder": "research/rank17_iv_surface_locus_certificate.py",
        },
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", type=Path, required=True)
    parser.add_argument("--small-prime", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build_certificate(arguments.production, arguments.small_prime)
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
        "cross_check": payload["independent_small_prime_cross_check"],
        "certificate_sha256": payload["certificate_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
