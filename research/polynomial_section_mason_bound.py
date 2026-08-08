#!/usr/bin/env python3
"""Exact certificate for the polynomial-degree bound and maximal passport.

Consider

    Q(v)^2 - v^2 L(v)^3 = R(v),

with deg(L)=2k, deg(Q)=3k+1, deg(R)<=3, Q(0)R(0)!=0.  Mason--
Stothers forces k<=2.  For k=2 and squarefree/coprime data, the associated
rational function Q^2/R is a genus-zero Belyi map with passport

    0:      2^7
    1:      3^4 2
    infinity: 11 1^3.

This script verifies every degree and Riemann--Hurwitz calculation exactly.
It records the theorem used by name; the mathematical proof is in the
accompanying document.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from research.degree24_polynomial_section_obstruction import canonical_payload


def mason_upper_radical_degree(k: int, common_degree: int = 0) -> int:
    """Upper bound for deg rad(A1 B1 C1) after removing the full gcd."""
    if k < 1:
        raise ValueError("k must be positive")
    if not 0 <= common_degree <= 3:
        raise ValueError("common_degree must lie between zero and three")
    # rad(Q^2/H) is supported on Q: at most 3k+1 roots.
    # rad(v^2 L^3/H) is supported on vL: at most 2k+1 roots.
    # rad(R/H) has degree at most 3-common_degree.
    return (3 * k + 1) + (2 * k + 1) + (3 - common_degree)


def maximal_term_degree(k: int, common_degree: int = 0) -> int:
    return 6 * k + 2 - common_degree


def verify_mason_inequality() -> None:
    for common_degree in range(4):
        for k in range(1, 12):
            radical = mason_upper_radical_degree(k, common_degree)
            maximum = maximal_term_degree(k, common_degree)
            permitted = maximum <= radical - 1
            if permitted != (k <= 2):
                raise AssertionError(
                    f"unexpected Mason outcome for k={k}, h={common_degree}"
                )


def cycle_ramification(cycles: list[int]) -> int:
    return sum(length - 1 for length in cycles)


def verify_passport() -> None:
    degree = 14
    zero_cycles = [2] * 7
    one_cycles = [3] * 4 + [2]
    infinity_cycles = [11, 1, 1, 1]

    for cycles in (zero_cycles, one_cycles, infinity_cycles):
        if sum(cycles) != degree:
            raise AssertionError("cycle partition has wrong degree")

    ramification = sum(
        cycle_ramification(cycles)
        for cycles in (zero_cycles, one_cycles, infinity_cycles)
    )
    if ramification != 2 * degree - 2:
        raise AssertionError("Riemann--Hurwitz genus-zero check failed")

    if cycle_ramification(zero_cycles) != 7:
        raise AssertionError("wrong ramification above zero")
    if cycle_ramification(one_cycles) != 9:
        raise AssertionError("wrong ramification above one")
    if cycle_ramification(infinity_cycles) != 10:
        raise AssertionError("wrong ramification above infinity")


def discriminant_coefficients() -> list[int]:
    """Return low-to-high coefficients of Disc(v^3-Sv^2+3v+1)."""
    return [-135, -54, 9, 4]


def multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def verify_cubic_degeneracies() -> None:
    # 4S^3+9S^2-54S-135=(S+3)^2(4S-15).
    expected = discriminant_coefficients()
    observed = multiply(multiply([3, 1], [3, 1]), [-15, 4])
    if observed != expected:
        raise AssertionError("cubic discriminant factorization failed")

    # Store polynomials low degree first and verify the two repeated-root cases.
    # S=-3: R=(v+1)^3.
    if multiply(multiply([1, 1], [1, 1]), [1, 1]) != [1, 3, 3, 1]:
        raise AssertionError("S=-3 factorization failed")
    # S=15/4: 4R=(v-2)^2(4v+1).
    if multiply(multiply([-2, 1], [-2, 1]), [1, 4]) != [4, -15, 12, 4]:
        raise AssertionError("S=15/4 factorization failed")


def verify_wronskian_degree() -> None:
    # In the squarefree/coprime k=2 case the Wronskian is divisible by
    # v*L^2*Q of degree 1+8+7=16, while rewriting it with R bounds its degree
    # by deg(Q^2)+deg(R)-1=14+3-1=16.  The quotient is therefore constant.
    divisor_degree = 1 + 2 * 4 + 7
    wronskian_upper_degree = 14 + 3 - 1
    if divisor_degree != 16 or wronskian_upper_degree != 16:
        raise AssertionError("maximal Wronskian degree calculation failed")


def build_certificate() -> dict:
    verify_mason_inequality()
    verify_passport()
    verify_cubic_degeneracies()
    verify_wronskian_degree()
    return {
        "certificate_id": "polynomial-section-mason-bound-v1",
        "claim_status": "proved",
        "degree_setup": {
            "deg_L": "2k",
            "deg_Q": "3k+1",
            "deg_R_upper_bound": 3,
            "identity": "Q(v)^2-v^2 L(v)^3=R(v)",
            "nonzero_at_zero": "Q(0)R(0)!=0",
        },
        "mason_stothers": {
            "after_full_gcd_H": {
                "deg_H": "h<=3",
                "maximum_term_degree": "6k+2-h",
                "radical_degree_upper_bound": "5k+5-h",
            },
            "inequality": "6k+2-h <= 5k+4-h",
            "result": "k<=2",
        },
        "maximal_case": {
            "degrees": {"L": 4, "Q": 7, "R": 3},
            "generic_assumptions": [
                "Q, vL, and R have disjoint supports",
                "Q, vL, and R are squarefree on their supports",
            ],
            "belyi_function": "f=Q^2/R",
            "passport": {
                "0": [2, 2, 2, 2, 2, 2, 2],
                "1": [3, 3, 3, 3, 2],
                "infinity": [11, 1, 1, 1],
            },
            "riemann_hurwitz": {
                "degree": 14,
                "ramification": 26,
                "required_for_genus_zero": 26,
            },
            "wronskian_identity": "2LQ+3vL'Q-2vLQ'=lambda, lambda in Q^*",
        },
        "target_cubic": {
            "R": "v^3-Sv^2+3v+1",
            "discriminant": "(S+3)^2(4S-15)",
            "degenerate_parameters": {
                "-3": "R=(v+1)^3",
                "15/4": "R=(v-2)^2(v+1/4)",
            },
        },
        "schema_version": 1,
        "scope_limitation": "The theorem bounds only this polynomial decomposition ansatz; rational functions with poles are not covered.",
        "theorem": "No polynomial identity in the stated degree pattern exists for k>=3; the sole unclassified higher level is (deg L,deg Q)=(4,7).",
        "verification": {
            "degree_inequalities": "exact integer arithmetic",
            "passport": "exact cycle partitions and Riemann--Hurwitz",
            "target_discriminant": "exact integer polynomial factorization",
            "wronskian_degree": "exact integer degree calculation",
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
    print("VERIFIED polynomial-section Mason bound and maximal passport")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
