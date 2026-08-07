#!/usr/bin/env sage-python
"""Rigorous ranks of the cyclic cubic bases u^3=c*t*(t-1)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, QQ, proof


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scales", default="1,2,4,9,12,18,27,36,48,108")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    proof.all(True)

    rows = []
    for token in args.scales.split(","):
        c = QQ(token.strip())
        # Omitting an explicit base-field positional argument keeps Sage in the
        # EllipticCurve_rational_field implementation rather than the generic
        # number-field implementation.
        curve = EllipticCurve([0, 0, 0, 0, 16 * c**4]).global_minimal_model()
        rank = curve.rank(proof=True)
        torsion = curve.torsion_subgroup()
        generators = curve.gens(proof=True) if rank > 0 else []
        rows.append(
            {
                "scale": str(c),
                "base_equation": f"u^3={c}*t*(t-1)",
                "mordell_model_a_invariants": [str(value) for value in curve.a_invariants()],
                "exact_rank": int(rank),
                "torsion_order": int(torsion.order()),
                "torsion_invariants": [int(value) for value in torsion.invariants()],
                "generators": [str(point) for point in generators],
            }
        )

    result = {
        "status": "pass",
        "transformation": "X=4*c*u, Y=4*c^2*(2*t-1), so Y^2=X^3+16*c^4",
        "rows": rows,
        "positive_rank_scales": [row["scale"] for row in rows if row["exact_rank"] > 0],
        "rank_zero_scales": [row["scale"] for row in rows if row["exact_rank"] == 0],
        "truth_status": "rigorous Sage ranks with proof.all(True) and rank(proof=True)",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "positive_rank_scales": result["positive_rank_scales"]}, sort_keys=True))


if __name__ == "__main__":
    main()
