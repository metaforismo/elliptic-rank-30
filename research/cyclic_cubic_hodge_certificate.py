#!/usr/bin/env python3
"""Exact arithmetic certificate for the C3 Hodge representation."""
from __future__ import annotations

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "cyclic_cubic_hodge_representation.json"

chi = 3
base_genus = 1
q = base_genus
p_g = chi + base_genus - 1
euler = 12 * chi
b2 = euler - 2 + 4 * q
h11 = b2 - 2 * p_g

assert (q, p_g, euler, b2, h11) == (1, 3, 36, 38, 32)

# Topological Lefschetz: fixed locus = three smooth elliptic fibres, Euler 0.
trace_h1 = -1
trace_h3 = -1
trace_h2 = -(1 - trace_h1 - trace_h3 + 1)
assert trace_h2 == -4

# Solve n0+2n=b2 and n0-n=trace_h2.
n = (b2 - trace_h2) // 3
n0 = trace_h2 + n
assert (n0, n) == (10, 14)
assert n0 + 2 * n == b2
assert n0 - n == trace_h2

# H20 multiplicities are (0,2,1), up to swapping the nontrivial characters.
h20 = (0, 2, 1)
h02 = (0, 1, 2)
h11_eigenspaces = tuple(
    total - holomorphic - antiholomorphic
    for total, holomorphic, antiholomorphic in zip((n0, n, n), h20, h02)
)
assert h11_eigenspaces == (10, 11, 11)

base_mw_rank = 8
nontrivial_capacity = h11_eigenspaces[1]
rank_ceiling = base_mw_rank + 2 * nontrivial_capacity
assert rank_ceiling == 30
picard_for_rank30 = rank_ceiling + 2
assert picard_for_rank30 == h11
transcendental_rank = b2 - h11
assert transcendental_rank == 2 * p_g == 6

certificate = {
    "status": "pass",
    "surface_invariants": {
        "chi": chi,
        "base_genus": base_genus,
        "q": q,
        "p_g": p_g,
        "euler_number": euler,
        "b2": b2,
        "h11": h11,
    },
    "deck_trace": {
        "fixed_locus_euler_number": 0,
        "trace_H1": trace_h1,
        "trace_H2": trace_h2,
        "trace_H3": trace_h3,
    },
    "H2_character_dimensions": {
        "trivial": n0,
        "zeta3": n,
        "zeta3_squared": n,
    },
    "H20_character_dimensions_up_to_conjugate_swap": {
        "trivial": h20[0],
        "zeta3": h20[1],
        "zeta3_squared": h20[2],
    },
    "H11_character_dimensions": {
        "trivial": h11_eigenspaces[0],
        "zeta3": h11_eigenspaces[1],
        "zeta3_squared": h11_eigenspaces[2],
    },
    "mordell_weil": {
        "invariant_rank": base_mw_rank,
        "nontrivial_character_capacity_each": nontrivial_capacity,
        "rank_ceiling": rank_ceiling,
        "rank30_equivalent_picard_number": picard_for_rank30,
        "transcendental_rank_at_picard_maximality": transcendental_rank,
    },
    "truth_status": "new intermediate theorem; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(
    json.dumps(
        {
            "status": certificate["status"],
            "H2": [n0, n, n],
            "H11": list(h11_eigenspaces),
            "rank_ceiling": rank_ceiling,
        },
        sort_keys=True,
    )
)
