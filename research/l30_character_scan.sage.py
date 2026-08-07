#!/usr/bin/env sage-python
"""Scan GAP character tables related to the automorphism group of the
30-dimensional 3-modular lattice for order-three elements with an
8-dimensional fixed space.

A cyclic cubic pullback of a split E8 rational elliptic surface requires an
order-three isometry of the putative rank-30 Mordell--Weil lattice whose fixed
lattice has rank 8. For a degree-30 character chi and an element g of order 3,
this fixed multiplicity is (chi(1)+chi(g)+chi(g^2))/3.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import ZZ, libgap

TABLE_NAMES = [
    "3_2.U4(3).(2^2)_133",
    "3_2.U4(3).2_3",
    "3_2.U4(3).2_1",
]


def sage_value(x):
    try:
        y = x.sage()
        if y in ZZ:
            return int(y)
        return str(y)
    except Exception:
        return str(x)


def scan_table(name):
    table = libgap.CharacterTable(name)
    irr = libgap.Irr(table)
    orders = [int(x) for x in libgap.OrdersClassRepresentatives(table)]
    names = [str(x) for x in libgap.ClassNames(table)]
    square_map = [int(x) - 1 for x in libgap.PowerMap(table, 2)]
    rows = []
    for character_index, chi in enumerate(irr, 1):
        degree = int(chi[0])
        if degree != 30:
            continue
        classes = []
        for i, order in enumerate(orders):
            if order != 3:
                continue
            value = chi[i]
            value_square = chi[square_map[i]]
            fixed = (chi[0] + value + value_square) / 3
            classes.append({
                "class_index": i + 1,
                "class_name": names[i],
                "character_value": sage_value(value),
                "square_class_index": square_map[i] + 1,
                "square_class_name": names[square_map[i]],
                "square_character_value": sage_value(value_square),
                "fixed_multiplicity": sage_value(fixed),
                "matches_cyclic_cubic_requirement": bool(fixed == 8),
            })
        rows.append({
            "character_index": character_index,
            "degree": degree,
            "order_three_classes": classes,
            "has_fixed_rank_8_class": any(
                c["matches_cyclic_cubic_requirement"] for c in classes
            ),
        })
    return {
        "table": name,
        "group_order": str(libgap.Size(table)),
        "class_count": len(orders),
        "degree_30_characters": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = {
        "status": "pass",
        "tables": [scan_table(name) for name in TABLE_NAMES],
        "truth_note": (
            "Character-table compatibility only. A fixed multiplicity of 8 does "
            "not by itself construct an integral lattice automorphism, an elliptic "
            "surface, or a rank-30 Mordell--Weil group."
        ),
    }
    out["matching_entries"] = [
        {
            "table": table["table"],
            "character_index": char["character_index"],
            "classes": [
                c for c in char["order_three_classes"]
                if c["matches_cyclic_cubic_requirement"]
            ],
        }
        for table in out["tables"]
        for char in table["degree_30_characters"]
        if char["has_fixed_rank_8_class"]
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": out["status"],
        "degree_30_character_count": sum(
            len(t["degree_30_characters"]) for t in out["tables"]
        ),
        "fixed_rank_8_match_count": len(out["matching_entries"]),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
