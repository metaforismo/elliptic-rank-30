#!/usr/bin/env sage-python
"""Classify the ADE root systems along a serialized rank-17 neighbor chain."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from sage.all import Genus, IntegralLattice, Matrix, QQ, ZZ, matrix, vector

ROOT_TYPES: dict[tuple[int, int], tuple[str, int]] = {}
for rank in range(1, 18):
    ROOT_TYPES[(rank, rank * (rank + 1))] = (f"A{rank}", rank + 1)
for rank in range(4, 18):
    ROOT_TYPES[(rank, 2 * rank * (rank - 1))] = (f"D{rank}", 4)
ROOT_TYPES[(6, 72)] = ("E6", 3)
ROOT_TYPES[(7, 126)] = ("E7", 2)
ROOT_TYPES[(8, 240)] = ("E8", 1)

FIBER_TYPES = {
    **{f"A{rank}": f"I{rank + 1}" for rank in range(2, 18)},
    "A1": "I2 or III",
    **{f"D{rank}": f"I{rank - 4}*" for rank in range(4, 18)},
    "E6": "IV*",
    "E7": "III*",
    "E8": "II*",
}


def serial_matrix(M: Matrix) -> list[list[int]]:
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


def union_find_components(roots: list, gram: Matrix) -> list[list[int]]:
    parent = list(range(len(roots)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    pairings = [root * gram for root in roots]
    for left in range(len(roots)):
        for right in range(left + 1, len(roots)):
            if pairings[left] * roots[right] != 0:
                union(left, right)

    groups: dict[int, list[int]] = {}
    for index in range(len(roots)):
        groups.setdefault(find(index), []).append(index)
    return sorted(groups.values(), key=lambda group: (len(group), group[0]), reverse=True)


def classify_entry(entry: dict, target_genus) -> dict:
    gram = Matrix(ZZ, entry["gram"])
    if gram.det() != 948 or Genus(gram) != target_genus:
        raise AssertionError(("entry outside target genus", entry["hash"], gram.det(), Genus(gram)))
    lattice = IntegralLattice(gram)
    roots = [vector(ZZ, row) for row in lattice.short_vectors(3)[2]]
    if len(roots) != entry["root_count"]:
        raise AssertionError(("root count mismatch", entry["hash"], len(roots), entry["root_count"]))

    components = []
    root_determinant = 1
    total_rank = 0
    if roots:
        for indices in union_find_components(roots, gram):
            vectors = matrix(ZZ, [list(roots[index]) for index in indices])
            rank = int(vectors.rank())
            count = len(indices)
            if (rank, count) not in ROOT_TYPES:
                raise AssertionError(("unclassified simply-laced component", rank, count, entry["hash"]))
            name, determinant = ROOT_TYPES[(rank, count)]
            total_rank += rank
            root_determinant *= determinant
            components.append({
                "type": name,
                "rank": rank,
                "root_count": count,
                "determinant": determinant,
                "kodaira_fiber_candidates": FIBER_TYPES[name],
            })

        root_module = matrix(ZZ, [list(root) for root in roots]).row_module(ZZ)
        root_basis = root_module.basis_matrix()
        saturation = root_module.saturation()
        saturation_index = int(root_module.index_in(saturation))
        root_gram = root_basis * gram * root_basis.transpose()
        if abs(int(root_gram.det())) != root_determinant:
            raise AssertionError(("root determinant mismatch", root_gram.det(), root_determinant, components))
    else:
        root_basis = matrix(ZZ, 0, 17)
        saturation_index = 1
        root_gram = matrix(ZZ, 0, 0)

    mw_rank = 17 - total_rank
    mw_regulator = QQ(948 * saturation_index * saturation_index) / QQ(root_determinant)
    return {
        "hash": entry["hash"],
        "source": entry.get("source"),
        "move": entry.get("move"),
        "gram": entry["gram"],
        "root_count": len(roots),
        "norm4_count": entry["norm4_count"],
        "root_system": " + ".join(component["type"] for component in components) if components else "0",
        "root_components": components,
        "root_rank": total_rank,
        "root_lattice_determinant": root_determinant,
        "root_lattice_primitive_closure_index": saturation_index,
        "torsion_order_predicted_by_shioda_tate": saturation_index,
        "mordell_weil_rank": mw_rank,
        "mordell_weil_regulator": str(mw_regulator),
        "rank_one_generator_height": str(mw_regulator) if mw_rank == 1 else None,
        "root_basis": serial_matrix(root_basis),
        "root_gram": serial_matrix(root_gram),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("chain", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    chain = json.loads(args.chain.read_text())
    period = Matrix(ZZ, [[-316, 0, 288], [0, 474, -15], [288, -15, -262]])
    target_genus = IntegralLattice(period).discriminant_group().genus((17, 0))
    records = [classify_entry(entry, target_genus) for entry in chain]

    seed = records[0]
    if seed["root_system"] != "A11 + A3 + A2":
        raise AssertionError(("unexpected transparent seed root system", seed["root_system"]))
    if seed["mordell_weil_rank"] != 1 or seed["rank_one_generator_height"] != "79/12":
        raise AssertionError(("unexpected seed Mordell-Weil data", seed))
    target = records[-1]
    if target["root_system"] != "0" or target["mordell_weil_rank"] != 17:
        raise AssertionError(("terminal lattice is not the rootless rank-17 target", target))

    payload = {
        "status": "completed",
        "claim_status": "exact ADE classification of every lattice in the frozen neighbor chain",
        "chain_sha256": hashlib.sha256(args.chain.read_bytes()).hexdigest(),
        "target_genus": str(target_genus),
        "records": records,
        "seed_fibration_data": {
            "essential_root_system": seed["root_system"],
            "reducible_fiber_root_lattice": seed["root_system"],
            "semistable_kodaira_configuration": ["I12", "I4", "I3"],
            "mordell_weil_rank": seed["mordell_weil_rank"],
            "torsion_order": seed["torsion_order_predicted_by_shioda_tate"],
            "generator_height": seed["rank_one_generator_height"],
            "remaining_euler_number_for_irreducible_singular_fibers": 5,
            "generic_semistable_configuration_candidate": "I12 + I4 + I3 + 5 I1",
        },
        "terminal_data": {
            "root_system": target["root_system"],
            "mordell_weil_rank": target["mordell_weil_rank"],
            "mordell_weil_regulator": target["mordell_weil_regulator"],
        },
        "truth_note": "The ADE and lattice computations are exact. Naming I2 rather than III for an isolated A1 would require the corresponding Weierstrass model; the transparent seed contains only A-type components of rank at least two, so its I12/I4/I3 labels are unambiguous.",
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "record_sha256": payload["record_sha256"],
        "root_systems": [record["root_system"] for record in records],
        "mw_ranks": [record["mordell_weil_rank"] for record in records],
        "seed_generator_height": seed["rank_one_generator_height"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
