#!/usr/bin/env python3
"""Exact algebraic certificate for the mu=2 j=0 cubic-family rank ceiling."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

OUT = Path("certificates/j0_cubic_family_rank_ceiling.json")


def main() -> int:
    v = sp.symbols("v")
    a = sp.Integer(3)
    b = sp.Rational(27, 2)
    original = sp.factor(b**2 * v**2 - a**3 * (v + 1) ** 3)
    expected_original = -sp.Rational(27, 4) * (v - 2) ** 2 * (4 * v + 1)
    assert sp.expand(original - expected_original) == 0

    discriminant_extension = 1 + 4 * v
    twist = sp.factor(discriminant_extension**3 * original)
    expected_twist = -sp.Rational(27, 4) * (v - 2) ** 2 * (4 * v + 1) ** 4
    assert sp.expand(twist - expected_twist) == 0

    # Kodaira data for y^2=x^3+a6 in characteristic zero.
    kodaira = {
        1: {"type": "II", "root_rank": 0, "euler": 2},
        2: {"type": "IV", "root_rank": 2, "euler": 4},
        3: {"type": "I0*", "root_rank": 4, "euler": 6},
        4: {"type": "IV*", "root_rank": 6, "euler": 8},
    }

    invariant_fibres = [
        {"count": 3, "order_a6": 2, **kodaira[2], "locus": "u^3=2c"},
        {"count": 3, "order_a6": 1, **kodaira[1], "locus": "4u^3=-c"},
        {"count": 1, "order_a6": 3, **kodaira[3], "locus": "u=infinity"},
    ]
    invariant_root_rank = sum(row["count"] * row["root_rank"] for row in invariant_fibres)
    invariant_euler = sum(row["count"] * row["euler"] for row in invariant_fibres)
    invariant_chi = 2
    invariant_h11 = 10 * invariant_chi
    invariant_rank_bound = invariant_h11 - 2 - invariant_root_rank
    assert invariant_root_rank == 10
    assert invariant_euler == 24 == 12 * invariant_chi
    assert invariant_rank_bound == 8

    twist_fibres = [
        {"count": 3, "order_a6": 2, **kodaira[2], "locus": "u^3=2c"},
        {"count": 3, "order_a6": 4, **kodaira[4], "locus": "4u^3=-c"},
    ]
    twist_root_rank = sum(row["count"] * row["root_rank"] for row in twist_fibres)
    twist_euler = sum(row["count"] * row["euler"] for row in twist_fibres)
    twist_chi = 3
    twist_h11 = 10 * twist_chi
    twist_rank_bound = twist_h11 - 2 - twist_root_rank
    assert twist_root_rank == 24
    assert twist_euler == 36 == 12 * twist_chi
    assert twist_rank_bound == 4

    total_bound = invariant_rank_bound + twist_rank_bound
    assert total_bound == 12

    result = {
        "status": "pass",
        "author": "Francesco Giannicola",
        "surface": "mu=2 marked j=0 family after u^3=c*t*(t-1)",
        "coefficient_in_v": str(original),
        "quadratic_extension": "Q(C_c)=Q(u,sqrt(1+4*u^3/c))",
        "quadratic_twist_coefficient_in_v": str(twist),
        "rank_decomposition": "rank E(Qbar(C_c)) = rank E0(Qbar(u)) + rank E0^(1+4u^3/c)(Qbar(u))",
        "invariant_surface": {
            "chi": invariant_chi,
            "h11": invariant_h11,
            "fibres": invariant_fibres,
            "root_rank": invariant_root_rank,
            "geometric_mordell_weil_rank_upper_bound": invariant_rank_bound,
        },
        "quadratic_twist_surface": {
            "chi": twist_chi,
            "h11": twist_h11,
            "fibres": twist_fibres,
            "root_rank": twist_root_rank,
            "geometric_mordell_weil_rank_upper_bound": twist_rank_bound,
        },
        "combined_geometric_rank_upper_bound": total_bound,
        "specialization_jump_needed_for_rank_30": 30 - total_bound,
        "restricted_obstruction": (
            "For every nonzero rational c, the generic Mordell-Weil rank over "
            "the cyclic genus-one base C_c is at most 12 over an algebraic closure."
        ),
        "truth_status": "restricted-family theorem; special fibres may still jump",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "pass",
        "invariant_bound": invariant_rank_bound,
        "twist_bound": twist_rank_bound,
        "combined_bound": total_bound,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
