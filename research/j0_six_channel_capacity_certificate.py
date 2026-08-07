#!/usr/bin/env python3
"""Exact six-character capacity certificate for the generic marked j=0 packet."""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

OUT = Path("certificates/j0_six_channel_capacity.json")


def kodaira(order: int) -> dict[str, int | str]:
    table = {
        0: {"type": "smooth", "root_rank": 0, "euler": 0},
        1: {"type": "II", "root_rank": 0, "euler": 2},
        2: {"type": "IV", "root_rank": 2, "euler": 4},
        3: {"type": "I0*", "root_rank": 4, "euler": 6},
        4: {"type": "IV*", "root_rank": 6, "euler": 8},
        5: {"type": "II*", "root_rank": 8, "euler": 10},
    }
    return table[order]


def channel(i: int, epsilon: int) -> dict:
    degree = 3 + 2 * i + 3 * epsilon
    chi = (degree + 5) // 6
    infinity_order = 6 * chi - degree
    finite = [
        {"locus": "three roots of P_mu(v)", "count": 3, "order_a6": 1, **kodaira(1)},
    ]
    if i:
        finite.append({"locus": "v=0", "count": 1, "order_a6": 2 * i, **kodaira(2 * i)})
    if epsilon:
        finite.append({"locus": "1+4v=0", "count": 1, "order_a6": 3, **kodaira(3)})
    if infinity_order:
        finite.append({"locus": "v=infinity", "count": 1, "order_a6": infinity_order, **kodaira(infinity_order)})
    root_rank = sum(row["count"] * row["root_rank"] for row in finite)
    euler = sum(row["count"] * row["euler"] for row in finite)
    h11 = 10 * chi
    rank_capacity = h11 - 2 - root_rank
    assert euler == 12 * chi
    return {
        "cubic_character": i,
        "quadratic_character": epsilon,
        "a6_multiplier": f"v^{2*i}*(1+4v)^{3*epsilon}",
        "degree_a6": degree,
        "chi": chi,
        "h11": h11,
        "infinity_order_a6": infinity_order,
        "fibres": finite,
        "root_rank": root_rank,
        "geometric_rank_capacity": rank_capacity,
        "picard_maximality_automatic": chi == 1,
    }


def main() -> int:
    mu, v = sp.symbols("mu v")
    a = mu**2 - mu + 1
    rmu = mu * (mu - 1)
    b = a**3 / rmu
    P = sp.factor(b**2 * v**2 - a**3 * (v + 1) ** 3)
    expected = sp.factor(
        -a**3
        * (v - mu * (mu - 1))
        * (mu**2 * v + mu - 1)
        * ((mu - 1) ** 2 * v - mu)
        / (mu**2 * (mu - 1) ** 2)
    )
    assert sp.simplify(P - expected) == 0

    roots = [
        sp.factor(mu * (mu - 1)),
        sp.factor((1 - mu) / mu**2),
        sp.factor(mu / (mu - 1) ** 2),
    ]
    assert sp.factor(roots[0] * roots[1] * roots[2]) == -1
    assert sp.factor(sum(roots[i] * roots[j] for i in range(3) for j in range(i + 1, 3))) == 3

    discriminant = sp.factor(
        sp.discriminant(sp.together(P * mu**2 * (mu - 1) ** 2 / (-a**3)), v)
    )
    expected_discriminant = sp.factor(a**6 * (a - 3) ** 2 * (a - 1) ** 2 * (4 * a - 3))
    assert sp.factor(discriminant - expected_discriminant) == 0
    collision_with_D = sp.factor(P.subs(v, -sp.Rational(1, 4)))
    assert sp.factor(collision_with_D) == sp.factor(
        a**3 * (a - 3) ** 2 * (4 * a - 3) / (64 * (a - 1) ** 2)
    )

    channels = [channel(i, epsilon) for epsilon in (0, 1) for i in (0, 1, 2)]
    capacities = [row["geometric_rank_capacity"] for row in channels]
    assert capacities == [4, 6, 4, 4, 6, 6]
    assert sum(capacities) == 30

    rational_channels = [row for row in channels if row["chi"] == 1]
    k3_channels = [row for row in channels if row["chi"] == 2]
    assert [row["geometric_rank_capacity"] for row in rational_channels] == [4, 6, 4]
    assert [row["geometric_rank_capacity"] for row in k3_channels] == [4, 6, 6]

    result = {
        "status": "pass",
        "author": "Francesco Giannicola",
        "family": "generic marked j=0 E8 family with cyclic cubic genus-one base",
        "excluded_markings": ["0", "1", "1/2", "2", "-1"],
        "P_mu_factorization": str(P),
        "three_rational_roots": [str(root) for root in roots],
        "P_mu_discriminant": str(discriminant),
        "P_mu_at_minus_one_quarter": str(collision_with_D),
        "character_decomposition": (
            "Over Qbar(v), the C2 x C3 packet is the direct sum of the six "
            "sextic-twist channels with a6 multipliers v^(2i)*(1+4v)^(3epsilon)."
        ),
        "channels": channels,
        "capacity_vector": capacities,
        "total_geometric_rank_capacity": sum(capacities),
        "automatic_rational_surface_rank_sum": sum(row["geometric_rank_capacity"] for row in rational_channels),
        "k3_capacity_vector": [row["geometric_rank_capacity"] for row in k3_channels],
        "rank_30_criterion": (
            "The three K3 channels must simultaneously have Picard number 20; "
            "the three rational elliptic channels are geometrically Picard maximal automatically."
        ),
        "truth_status": "new intermediate theorem; rational descent and explicit sections remain open",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "pass",
        "capacity_vector": capacities,
        "total": sum(capacities),
        "k3_capacities": result["k3_capacity_vector"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
