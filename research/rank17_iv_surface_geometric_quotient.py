#!/usr/bin/env python3
"""Quotient the finite-field IV surface tuples by the jet-sign involution.

The transformation

    (e0,p0,p1,p2,p3,r,s)
      -> (-e0,-p0,-p1,-p2,-p3,r,s)

fixes c4, c6, the discriminant, the split tangent targets, and e0*P.  Hence it
changes only an auxiliary square/cube-jet presentation, not the Weierstrass
surface or the adapted section coordinate W=e0*P*D^2+t*U.

This script verifies every orbit exactly and keeps the e0=1 representative.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def parameter_tuple(candidate: dict[str, object]) -> tuple[int, ...]:
    params = candidate["parameters"]
    return tuple(
        int(params[name])
        for name in ("e0", "p0", "p1", "p2", "p3", "r", "s")
    )


def partner_tuple(values: tuple[int, ...], prime: int) -> tuple[int, ...]:
    e0, p0, p1, p2, p3, r, s = values
    return (
        (-e0) % prime,
        (-p0) % prime,
        (-p1) % prime,
        (-p2) % prime,
        (-p3) % prime,
        r,
        s,
    )


def invariant_projection(candidate: dict[str, object]) -> dict[str, object]:
    return {
        "c4_coefficients_ascending": candidate["c4_coefficients_ascending"],
        "c6_coefficients_ascending": candidate["c6_coefficients_ascending"],
        "discriminant_coefficients_ascending": candidate[
            "discriminant_coefficients_ascending"
        ],
        "exact_fibre_checks": candidate["exact_fibre_checks"],
        "split_tangent_checks": candidate["split_tangent_checks"],
    }


def build(input_path: Path) -> dict[str, object]:
    source = json.loads(input_path.read_text(encoding="utf-8"))
    prime = int(source["prime"])
    candidates = source["candidates"]
    by_tuple = {parameter_tuple(candidate): candidate for candidate in candidates}
    if len(by_tuple) != len(candidates):
        raise AssertionError("duplicate raw surface tuples")

    visited: set[tuple[int, ...]] = set()
    orbits = []
    representatives = []
    for values in sorted(by_tuple):
        if values in visited:
            continue
        partner_values = partner_tuple(values, prime)
        if partner_values == values:
            raise AssertionError(("unexpected fixed point", values))
        if partner_values not in by_tuple:
            raise AssertionError(("missing involution partner", values, partner_values))
        left = by_tuple[values]
        right = by_tuple[partner_values]
        if invariant_projection(left) != invariant_projection(right):
            raise AssertionError(("involution changed the geometric surface", values))
        visited.add(values)
        visited.add(partner_values)

        representative_values = values if values[0] == 1 else partner_values
        representative = by_tuple[representative_values]
        if representative_values[0] != 1:
            raise AssertionError("the orbit has no e0=1 representative")
        orbit = {
            "raw_parameter_tuples": [list(values), list(partner_values)],
            "representative_parameter_tuple": list(representative_values),
            "representative_record_sha256": representative["record_sha256"],
            "surface_invariant_sha256": hashlib.sha256(
                json.dumps(
                    invariant_projection(representative),
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
        }
        orbits.append(orbit)
        representatives.append(representative)

    if len(visited) != len(candidates):
        raise AssertionError("the involution orbits do not cover every candidate")
    if 2 * len(orbits) != len(candidates):
        raise AssertionError("an orbit has size different from two")

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": f"rank17_iv_surface_geometric_quotient_f{prime}",
        "truth_status": (
            f"EXACT quotient of the normalized F_{prime} surface tuples by the "
            "auxiliary jet-sign involution; no characteristic-zero or rank-30 conclusion"
        ),
        "prime": prime,
        "source_certificate_sha256": source["certificate_sha256"],
        "raw_candidate_count": len(candidates),
        "geometric_surface_count": len(orbits),
        "group_action": {
            "order": 2,
            "formula": (
                "(e0,p0,p1,p2,p3,r,s)->"
                "(-e0,-p0,-p1,-p2,-p3,r,s)"
            ),
            "fixed_geometric_data": [
                "c4",
                "c6",
                "discriminant",
                "split tangent targets",
                "e0*P",
            ],
            "every_orbit_has_size_two": True,
        },
        "orbits": orbits,
        "representatives": representatives,
        "section_sign_quotient_note": (
            "For a fixed surface, negating a section simultaneously negates the "
            "I4 tangent root, the IV tangent root, and q.  Thus the four local "
            "root choices form two further sign orbits; the first probe retains "
            "all four as an independent symmetry check."
        ),
        "limitations": [
            "This quotient removes only the explicit auxiliary jet-sign redundancy.",
            "It does not decide whether two different invariant records are isomorphic by another base or Weierstrass automorphism.",
            "Finite-field surfaces need not lift to characteristic zero.",
        ],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build(arguments.input)
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
        "raw_candidate_count": payload["raw_candidate_count"],
        "geometric_surface_count": payload["geometric_surface_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
