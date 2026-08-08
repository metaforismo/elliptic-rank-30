#!/usr/bin/env python3
"""Exact arithmetic certificate for the one-function visibility barrier.

Let E/k be a short Weierstrass curve y^2=R(x), deg R=3, with identity O.
Every f in L(N O) is represented as A(x)+y B(x), and the pole orders of
x and y at O are 2 and 3.  Because 2*deg(A) is even and 2*deg(B)+3 is odd,
the leading pole terms cannot cancel, so

    ord_O^-(A+yB) = max(2 deg A, 2 deg B+3).

If all N zeros of f are rational, the principal divisor gives one exact
relation among the corresponding points.  Hence one function with N=30
visible zeros can span rank at most 29; N=31 is the first compatible pole
order for a rank-30 span.

The script certifies the exact degree/pole-order arithmetic and the
minimal norm-form shape

    A(x)^2 - R(x) B(x)^2,

with deg A <= 15, deg B = 14, deg R = 3, whose generic degree is 31.
It does not claim that such a norm polynomial splits over Q or that the
remaining thirty point classes are independent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Optional


def pole_order(deg_a: Optional[int], deg_b: Optional[int]) -> int:
    """Pole order at O of A(x)+yB(x); None denotes the zero polynomial."""
    orders: list[int] = []
    if deg_a is not None:
        if deg_a < 0:
            raise ValueError("deg_a must be nonnegative or None")
        orders.append(2 * deg_a)
    if deg_b is not None:
        if deg_b < 0:
            raise ValueError("deg_b must be nonnegative or None")
        orders.append(2 * deg_b + 3)
    if not orders:
        raise ValueError("the zero function has no exact pole order")
    # The A and yB orders have opposite parity, so equality/cancellation is
    # impossible when both terms are nonzero.
    return max(orders)


def norm_degree(
    deg_a: Optional[int], deg_b: Optional[int], deg_r: int = 3
) -> int:
    """Generic degree of A^2-RB^2 when the top terms do not cancel."""
    degrees: list[int] = []
    if deg_a is not None:
        degrees.append(2 * deg_a)
    if deg_b is not None:
        degrees.append(deg_r + 2 * deg_b)
    if not degrees:
        raise ValueError("zero norm expression")
    return max(degrees)


def l_space_basis_pole_orders(n: int) -> list[int]:
    """Pole orders of the standard basis of L(nO) for n>=0."""
    if n < 0:
        return []
    orders = [0]
    orders.extend(range(2, n + 1, 2))
    orders.extend(range(3, n + 1, 2))
    return sorted(orders)


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compute_certificate() -> dict[str, object]:
    # Riemann--Roch on a genus-one curve gives ell(nO)=n for n>=1.
    dimensions: dict[str, int] = {}
    basis_orders: dict[str, list[int]] = {}
    for n in range(1, 32):
        orders = l_space_basis_pole_orders(n)
        if len(orders) != n:
            raise AssertionError(
                f"L({n}O) basis has {len(orders)} elements, expected {n}"
            )
        if orders[0] != 0 or orders[-1] > n:
            raise AssertionError(f"invalid basis pole orders for n={n}")
        dimensions[str(n)] = len(orders)
        basis_orders[str(n)] = orders

    # Verify the near-square mechanism for every nontrivial degree up to 15.
    near_square_table: list[dict[str, int]] = []
    for degree_q in range(2, 16):
        order = pole_order(degree_q, None)
        degree_of_q2_minus_r = max(2 * degree_q, 3)
        if order != 2 * degree_q:
            raise AssertionError("near-square pole-order formula failed")
        if degree_of_q2_minus_r != order:
            raise AssertionError("near-square norm degree mismatch")
        near_square_table.append(
            {
                "degree_Q": degree_q,
                "pole_order": order,
                "visible_zero_count_if_split_simple": order,
                "rank_upper_bound_from_visible_points": order - 1,
            }
        )

    if near_square_table[-1] != {
        "degree_Q": 15,
        "pole_order": 30,
        "visible_zero_count_if_split_simple": 30,
        "rank_upper_bound_from_visible_points": 29,
    }:
        raise AssertionError("degree-15 near-square frontier mismatch")

    # Enumerate all degree patterns that have exact pole order 30 or 31.
    patterns: dict[str, list[dict[str, Optional[int]]]] = {"30": [], "31": []}
    for deg_a in [None, *range(0, 17)]:
        for deg_b in [None, *range(0, 16)]:
            if deg_a is None and deg_b is None:
                continue
            order = pole_order(deg_a, deg_b)
            if order in (30, 31):
                patterns[str(order)].append(
                    {
                        "deg_A": deg_a,
                        "deg_B": deg_b,
                        "generic_norm_degree": norm_degree(deg_a, deg_b),
                    }
                )

    # Exact pole order 31 must come from yB with deg B=14; A can have deg<=15.
    expected_31 = [
        {
            "deg_A": deg_a,
            "deg_B": 14,
            "generic_norm_degree": 31,
        }
        for deg_a in [None, *range(0, 16)]
    ]
    if patterns["31"] != expected_31:
        raise AssertionError(
            f"unexpected order-31 patterns: {patterns['31']}"
        )

    # Exact pole order 30 comes from deg A=15 and deg B<=13 (or B=0).
    expected_30 = [
        {
            "deg_A": 15,
            "deg_B": deg_b,
            "generic_norm_degree": 30,
        }
        for deg_b in [None, *range(0, 14)]
    ]
    if patterns["30"] != expected_30:
        raise AssertionError(
            f"unexpected order-30 patterns: {patterns['30']}"
        )

    payload: dict[str, object] = {
        "schema_version": 1,
        "theorem": {
            "name": "one-function visibility barrier",
            "statement": (
                "If a nonzero rational function f on E has pole divisor N*O "
                "and N rational zeros counted with multiplicity, then their "
                "sum is O in E(k). Therefore those N visible points span rank "
                "at most N-1."
            ),
            "proof_object": (
                "(f)=sum_i(P_i)-N(O) is principal, and Pic^0(E) is E."
            ),
            "conditional_assumptions": [],
        },
        "riemann_roch_basis": {
            "curve_model": "y^2=R(x), deg R=3, char(k) not 2 or 3",
            "pole_orders": {"x": 2, "y": 3},
            "orders_for_L_30O": basis_orders["30"],
            "orders_for_L_31O": basis_orders["31"],
            "dimensions_1_through_31": dimensions,
            "representation": "f=A(x)+y*B(x)",
            "exact_pole_order": "max(2*deg(A),2*deg(B)+3)",
            "no_leading_cancellation_reason": (
                "the A and yB pole orders have opposite parity"
            ),
        },
        "near_square_frontier": {
            "function": "f=y-Q(x)",
            "identity": "Norm(f)=Q(x)^2-R(x)",
            "table_degrees_2_through_15": near_square_table,
            "degree_15_conclusion": (
                "Thirty split simple rational zeros carry the forced relation "
                "P1+...+P30=O, so their span has rank at most 29."
            ),
        },
        "minimal_single_function_rank30_target": {
            "minimal_pole_order": 31,
            "degree_patterns": patterns["31"],
            "canonical_choice": {
                "deg_A_max": 15,
                "deg_B": 14,
                "deg_R": 3,
                "norm_form": "A(x)^2-R(x)*B(x)^2",
                "generic_norm_degree": 31,
            },
            "sufficient_point_recovery_conditions": [
                "R has nonzero discriminant",
                "A^2-R*B^2 has 31 distinct rational roots r_i",
                "B(r_i) is nonzero for every i",
                "P_i=(r_i,-A(r_i)/B(r_i))",
            ],
            "forced_relation": "P1+...+P31=O",
            "maximum_possible_span_after_forced_relation": 30,
            "truth_note": (
                "The norm identity and splitting conditions are only a "
                "construction target. They do not prove that the remaining "
                "thirty classes are independent."
            ),
        },
        "order_pattern_audit": patterns,
        "implementation": {
            "language": "Python standard library",
            "script": "research/riemann_roch_visibility_certificate.py",
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the freshly recomputed JSON certificate",
    )
    parser.add_argument(
        "--compare",
        type=Path,
        help="compare recomputation with a committed JSON certificate",
    )
    arguments = parser.parse_args()

    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != certificate:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
        print(f"matched {arguments.compare}")

    print("L(30O) dimension: 30; L(31O) dimension: 31")
    print("degree-15 near-square visible span: at most 29")
    print("minimal one-function rank-30-compatible pole order: 31")
    print("canonical norm target: deg A<=15, deg B=14, deg Norm=31")
    print(f"certificate sha256: {certificate['certificate_sha256']}")


if __name__ == "__main__":
    main()
