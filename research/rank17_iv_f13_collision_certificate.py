#!/usr/bin/env sage -python
"""Certify the residual-pole collision of the unique split F13 IV surface."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.version import version as sage_version


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build(path: Path) -> dict[str, object]:
    source = json.loads(path.read_text(encoding="utf-8"))
    if int(source["prime"]) != 13:
        raise ValueError("the collision certificate is specific to F13")
    candidates = [
        candidate for candidate in source["candidates"]
        if int(candidate["parameters"]["e0"]) == 1
    ]
    if len(candidates) != 1:
        raise AssertionError(("unexpected canonical F13 candidate count", len(candidates)))
    candidate = candidates[0]
    field = GF(13)
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    c4 = ring(candidate["c4_coefficients_ascending"])
    c6 = ring(candidate["c6_coefficients_ascending"])
    A, rem4 = c4.quo_rem((t-1)**2)
    B, rem6 = c6.quo_rem((t-1)**2)
    if rem4 or rem6:
        raise AssertionError("the F13 invariants lost the IV factors")
    C = 2*A*B + 3*(t-1)*A.derivative()*B - 2*(t-1)*A*B.derivative()
    support = [index for index in range(C.degree()+1) if C[index]]
    if support != [3, 4]:
        raise AssertionError(("unexpected Wronskian support", support))
    extra_point = -C[3]/C[4]
    H = (t-1)**2*A**3 - B**2
    R, remainder = H.quo_rem(t**4)
    if remainder or R.degree() != 4:
        raise AssertionError("the residual quartic reconstruction failed")
    gcd = R.gcd(R.derivative()).monic()
    factorization = R.factor()
    if extra_point != field(8):
        raise AssertionError(("unexpected extra ramification point", extra_point))
    if R(extra_point) != 0:
        raise AssertionError("the extra point does not collide with a residual pole")
    if gcd != t + 5:
        raise AssertionError(("unexpected repeated residual factor", gcd))
    expected = (-3) * (t+5)**2 * (t**2 + 2*t + 6)
    if R != expected:
        raise AssertionError(("unexpected residual factorization", R.factor()))

    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_f13_extra_branch_pole_collision",
        "truth_status": (
            "CERTIFIED degeneration of the unique normalized split F13 surface; "
            "this is not a characteristic-zero obstruction or rank-30 conclusion"
        ),
        "sage_version": str(sage_version),
        "prime": 13,
        "source_surface_certificate_sha256": source["certificate_sha256"],
        "source_candidate_sha256": candidate["record_sha256"],
        "parameters": candidate["parameters"],
        "wronskian_C_coefficients_ascending": [int(value) for value in C.list()],
        "wronskian_support": support,
        "extra_ramification_point": int(extra_point),
        "residual_quartic_coefficients_ascending": [int(value) for value in R.list()],
        "residual_quartic_factorization": "-3*(t+5)^2*(t^2+2*t+6)",
        "gcd_residual_derivative": [int(value) for value in gcd.list()],
        "extra_point_is_repeated_residual_pole": True,
        "geometric_consequence": (
            "The sole split F13 point lies on the boundary where the fourth "
            "simple branch point collides with a pole. It is not in the open "
            "I12+I4+IV+4I1 locus."
        ),
        "limitations": [
            "A Q-surface may have bad or nonintegral reduction at 13.",
            "The result does not exclude the additive-IV locus over Q.",
            "No section is imposed and the global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
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
        "extra_ramification_point": payload["extra_ramification_point"],
        "residual_quartic_factorization": payload["residual_quartic_factorization"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
