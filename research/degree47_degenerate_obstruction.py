#!/usr/bin/env python3
"""Exact degree bookkeeping for the two degenerate maximal remainders.

For the target R(v)=v^3-Sv^2+3v+1, the repeated-root cases are

    S=-3:   R=(v+1)^3,
    S=15/4: R=(v-2)^2(v+1/4).

Local valuations plus Mason--Stothers rule out a degree-(4,7) identity in
both cases.  This is an intermediate obstruction, not a rank-30 certificate.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from research.degree24_polynomial_section_obstruction import canonical_payload


def verify_triple_root_case() -> None:
    # Branch 1: Q and L are units at v=-1.  Radical support degrees are at
    # most 7 (Q), 5 (vL), and 1 (R); Mason gives 14 <= 12, impossible.
    maximum = 14
    radical = 7 + 5 + 1
    if maximum <= radical - 1:
        raise AssertionError("the unit branch was not obstructed")

    # Branch 2: ord_{-1}(L)=1 and ord_{-1}(Q)=q>=2.  After removing R as
    # the full gcd, both large terms have degree 11.  The radical has at most
    # (1+7-q) roots from Q^2/R and 4 from v^2(L/(v+1))^3.
    for q in range(2, 8):
        radical_after_gcd = (1 + 7 - q) + 4
        if 11 <= radical_after_gcd - 1:
            raise AssertionError(f"triple-root branch survived for q={q}")


def verify_double_root_case() -> None:
    # Branch 1: Q and L are units at v=2.  R has two distinct roots.
    maximum = 14
    radical = 7 + 5 + 2
    if maximum <= radical - 1:
        raise AssertionError("the double-root unit branch was not obstructed")

    # Branch 2: ord_2(Q)=1 and ord_2(L)=l>=1.  Removing (v-2)^2 leaves
    # degree-12 large terms.  Radical support: at most 6 roots from Q/(v-2),
    # 2+(4-l) roots from v, v-2, and the remaining L factor, and one root
    # from v+1/4.
    for l in range(1, 5):
        radical_after_gcd = 6 + (2 + 4 - l) + 1
        if 12 <= radical_after_gcd - 1:
            raise AssertionError(f"double-root branch survived for l={l}")


def verify_local_order_classification() -> None:
    # Enumerate possible nonnegative local orders within the polynomial
    # degree ranges and confirm the only ways to obtain remainder order 3 or
    # 2, allowing the order-zero cancellation branch.
    triple_positive = []
    double_positive = []
    for q in range(8):
        for l in range(5):
            if q == l == 0:
                continue
            left = 2 * q
            right = 3 * l
            if left != right and min(left, right) == 3:
                triple_positive.append((q, l))
            if left != right and min(left, right) == 2:
                double_positive.append((q, l))
    if triple_positive != [(q, 1) for q in range(2, 8)]:
        raise AssertionError(f"unexpected triple-root local branches: {triple_positive}")
    if double_positive != [(1, l) for l in range(1, 5)]:
        raise AssertionError(f"unexpected double-root local branches: {double_positive}")


def build_certificate() -> dict:
    verify_local_order_classification()
    verify_triple_root_case()
    verify_double_root_case()
    return {
        "certificate_id": "degree47-degenerate-remainder-obstruction-v1",
        "claim_status": "proved",
        "degree_setup": {
            "L": 4,
            "Q": 7,
            "identity": "Q(v)^2-v^2L(v)^3=R(v)",
        },
        "double_root_case": {
            "parameter": "S=15/4",
            "remainder": "R=(v-2)^2(v+1/4)",
            "unit_branch": "Mason: 14 <= 7+5+2-1=13, contradiction",
            "common_branch": {
                "local_orders": "ord_2(Q)=1, ord_2(L)=l>=1",
                "after_gcd_degree": 12,
                "radical_upper_bound": "13-l",
                "mason_contradiction": "12 <= 12-l, impossible for l>=1",
            },
            "result": "no degree-(4,7) solution",
        },
        "schema_version": 1,
        "scope_limitation": "Rules out only the two repeated-root target cubics at the maximal polynomial degree.",
        "theorem": "Neither repeated-root value S=-3 nor S=15/4 admits a degree-(4,7) polynomial identity.",
        "triple_root_case": {
            "parameter": "S=-3",
            "remainder": "R=(v+1)^3",
            "unit_branch": "Mason: 14 <= 7+5+1-1=12, contradiction",
            "common_branch": {
                "local_orders": "ord_{-1}(L)=1, ord_{-1}(Q)=q>=2",
                "after_gcd_degree": 11,
                "radical_upper_bound": "12-q",
                "mason_contradiction": "11 <= 11-q, impossible for q>=2",
            },
            "result": "no degree-(4,7) solution",
        },
        "verification": {
            "local_order_branches": "complete finite enumeration of possible local orders",
            "mason_degree_bounds": "exact integer arithmetic",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    payload = canonical_payload(build_certificate())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check and args.check.read_text(encoding="utf-8") != payload:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print("VERIFIED both degenerate degree-(4,7) obstructions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
