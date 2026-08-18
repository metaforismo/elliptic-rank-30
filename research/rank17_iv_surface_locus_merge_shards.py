#!/usr/bin/env python3
"""Merge a complete disjoint p0-sharded additive-IV surface census.

Each input is a deterministic JSON output of
rank17_iv_surface_locus_all_shard.c.  The intervals must be disjoint,
contiguous, and cover [1,p).  Candidate tuples must lie in their declared
interval, be unique globally, and their split-category counts must replay.
The result has exactly the schema consumed by
rank17_iv_surface_locus_all_certificate.py.

This is repository plumbing for an exact finite-field exhaustion; it makes no
characteristic-zero, section, or rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

CATEGORY_FROM_FLAGS = {
    (0, 0): "neither",
    (0, 1): "iv_only",
    (1, 0): "i4_only",
    (1, 1): "both",
}


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "prime",
        "normalization",
        "p0_start_inclusive",
        "p0_end_exclusive",
        "visited_parameter_tuples",
        "candidate_count",
        "split_category_counts",
        "tuples",
    }
    missing = required - set(data)
    if missing:
        raise ValueError((path, sorted(missing)))
    return data


def build(paths: list[Path]) -> dict[str, object]:
    if not paths:
        raise ValueError("at least one shard is required")
    loaded = [(path, load(path)) for path in paths]
    primes = {int(data["prime"]) for _path, data in loaded}
    if len(primes) != 1:
        raise ValueError(("shard primes differ", sorted(primes)))
    prime = primes.pop()
    normalizations = {
        str(data["normalization"]) for _path, data in loaded
    }
    if normalizations != {"e0=1 jet-sign representative"}:
        raise ValueError(("unexpected normalizations", sorted(normalizations)))

    loaded.sort(key=lambda item: int(item[1]["p0_start_inclusive"]))
    expected_start = 1
    tuples: list[tuple[int, ...]] = []
    shard_records = []
    category_totals = {
        "neither": 0,
        "iv_only": 0,
        "i4_only": 0,
        "both": 0,
    }
    total_visited = 0

    for path, data in loaded:
        start = int(data["p0_start_inclusive"])
        end = int(data["p0_end_exclusive"])
        if start != expected_start:
            raise AssertionError((
                "p0 shards are not a contiguous partition",
                expected_start,
                start,
                path,
            ))
        if not (1 <= start < end <= prime):
            raise AssertionError(("invalid p0 interval", path, start, end))
        expected_start = end

        expected_visited = (end - start) * prime**5
        visited = int(data["visited_parameter_tuples"])
        if visited != expected_visited:
            raise AssertionError((
                "shard visited-count mismatch",
                path,
                visited,
                expected_visited,
            ))
        total_visited += visited

        local_tuples = []
        local_categories = {
            "neither": 0,
            "iv_only": 0,
            "i4_only": 0,
            "both": 0,
        }
        for raw in data["tuples"]:
            value = tuple(int(entry) for entry in raw)
            if len(value) != 8:
                raise AssertionError(("unexpected tuple length", path, value))
            p0, p1, p2, p3, r, s, split_i4, split_iv = value
            if not start <= p0 < end:
                raise AssertionError(("tuple escapes shard interval", path, value))
            if any(not 0 <= entry < prime for entry in (p0,p1,p2,p3,r,s)):
                raise AssertionError(("tuple coordinate outside F_p", path, value))
            if split_i4 not in (0, 1) or split_iv not in (0, 1):
                raise AssertionError(("invalid split flags", path, value))
            category = CATEGORY_FROM_FLAGS[(split_i4, split_iv)]
            local_categories[category] += 1
            local_tuples.append(value)

        if len(local_tuples) != int(data["candidate_count"]):
            raise AssertionError(("shard candidate count mismatch", path))
        if local_categories != {
            key: int(data["split_category_counts"][key])
            for key in local_categories
        }:
            raise AssertionError(("shard category replay mismatch", path))
        if len(set(local_tuples)) != len(local_tuples):
            raise AssertionError(("duplicate tuple inside shard", path))

        tuples.extend(local_tuples)
        for key in category_totals:
            category_totals[key] += local_categories[key]
        shard_records.append({
            "path": str(path),
            "p0_start_inclusive": start,
            "p0_end_exclusive": end,
            "visited_parameter_tuples": visited,
            "candidate_count": len(local_tuples),
            "split_category_counts": local_categories,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })

    if expected_start != prime:
        raise AssertionError((
            "p0 shards do not cover the full normalized range",
            expected_start,
            prime,
        ))
    if total_visited != (prime - 1) * prime**5:
        raise AssertionError((
            "global visited-count mismatch",
            total_visited,
            (prime - 1) * prime**5,
        ))
    tuples.sort()
    if len(set(tuples)) != len(tuples):
        raise AssertionError("duplicate candidate across shards")

    production = {
        "prime": prime,
        "normalization": "e0=1 jet-sign representative",
        "visited_parameter_tuples": total_visited,
        "candidate_count": len(tuples),
        "split_category_counts": category_totals,
        "tuples": [list(value) for value in tuples],
    }
    payload = {
        "schema_version": 1,
        "truth_status": (
            f"EXACT merge of a complete disjoint F_{prime} p0-shard "
            "exhaustion; no characteristic-zero, section, or rank-30 conclusion"
        ),
        "prime": prime,
        "partition": shard_records,
        "production": production,
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--production-output", type=Path, required=True)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build(arguments.shard)
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != payload:
            raise AssertionError(f"merge certificate mismatch: {arguments.compare}")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    arguments.production_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.production_output.write_text(
        json.dumps(payload["production"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "prime": payload["prime"],
        "shard_count": len(payload["partition"]),
        "visited_parameter_tuples": payload["production"]["visited_parameter_tuples"],
        "candidate_count": payload["production"]["candidate_count"],
        "split_category_counts": payload["production"]["split_category_counts"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
