#!/usr/bin/env sage-python
"""Fast exact E8(3)-embedding search in the extremal rank-30 3-modular lattice.

This imports the catalogue Gram matrix and certificate helpers from
``l30_e8_embedding.sage.py`` but replaces both bottlenecks:

* PARI ``qfminim`` enumerates the 13,320 minimal vectors;
* inner-product relation sets are represented as Python integer bitsets, and
  the long E8 arm is built before the two interchangeable short arms.

A positive result is a fully exact lattice theorem.  A negative result is only
complete when every minimal vector has been tried as the trivalent root.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

from sage.all import IntegralLattice, Matrix, ZZ, pari, vector

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "l30_slow", ROOT / "l30_e8_embedding.sage.py"
)
SLOW = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SLOW)


def enumerate_norm_six(G):
    """Return one integral vector from each pair +/-v with vGv^t=6."""
    pg = pari(G)
    try:
        result = pg.qfminim(ZZ(6), ZZ(7000), 0)
    except TypeError:
        result = pari.qfminim(pg, ZZ(6), ZZ(7000), 0)
    total = int(result[0])
    raw = Matrix(ZZ, result[2].sage())
    if raw.nrows() == G.nrows():
        rows = raw.transpose()
    elif raw.ncols() == G.nrows():
        rows = raw
    else:
        raise RuntimeError(f"unexpected qfminim matrix shape {raw.dimensions()}")
    vectors = [vector(ZZ, row) for row in rows.rows()]
    if total != 13320 or len(vectors) != 6660:
        raise AssertionError(
            f"expected 13320 vectors / 6660 sign pairs, got {total} / {len(vectors)}"
        )
    if any((v * G).dot_product(v) != 6 for v in vectors):
        raise AssertionError("qfminim returned a vector of norm different from 6")
    return vectors, total


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


class Relations:
    def __init__(self, vectors, G):
        self.vectors = vectors
        self.V = Matrix(ZZ, vectors)
        self.GV = self.V * G
        self.VT = self.V.transpose()
        self.cache = {}
        self.all_mask = (1 << len(vectors)) - 1

    def get(self, index):
        cached = self.cache.get(index)
        if cached is not None:
            return cached
        values = (self.GV.row(index) * self.VT).row(0)
        neighbor3 = 0
        orthogonal = 0
        for j, value in enumerate(values):
            value = int(value)
            if abs(value) == 3:
                neighbor3 |= 1 << j
            elif value == 0:
                orthogonal |= 1 << j
        cached = (neighbor3, orthogonal)
        self.cache[index] = cached
        return cached


def search_embedding(vectors, G, central_limit):
    rel = Relations(vectors, G)
    checked = 0
    used_central = 0
    limit = min(central_limit, len(vectors))

    for central in range(limit):
        used_central += 1
        n3_c, z_c = rel.get(central)

        # Build the length-five arm 2-3-4-5-6-7 first.  Each new vertex is
        # adjacent to the previous one and orthogonal to every earlier
        # non-neighbour.
        for node3 in bits(n3_c):
            checked += 1
            n3_3, z_3 = rel.get(node3)
            mask4 = n3_3 & z_c
            for node4 in bits(mask4):
                checked += 1
                n3_4, z_4 = rel.get(node4)
                mask5 = n3_4 & z_c & z_3
                for node5 in bits(mask5):
                    checked += 1
                    n3_5, z_5 = rel.get(node5)
                    mask6 = n3_5 & z_c & z_3 & z_4
                    for node6 in bits(mask6):
                        checked += 1
                        n3_6, z_6 = rel.get(node6)
                        mask7 = n3_6 & z_c & z_3 & z_4 & z_5
                        for node7 in bits(mask7):
                            checked += 1
                            _, z_7 = rel.get(node7)

                            # The two short arms are both adjacent to the
                            # central root and orthogonal to the complete long
                            # arm.  They must also be mutually orthogonal.
                            arms = n3_c & z_3 & z_4 & z_5 & z_6 & z_7
                            for node0 in bits(arms):
                                _, z_0 = rel.get(node0)
                                later = arms & z_0
                                # Remove the symmetry interchanging nodes 0,1.
                                later &= ~((1 << (node0 + 1)) - 1)
                                for node1 in bits(later):
                                    checked += 1
                                    return (
                                        {
                                            0: node0,
                                            1: node1,
                                            2: central,
                                            3: node3,
                                            4: node4,
                                            5: node5,
                                            6: node6,
                                            7: node7,
                                        },
                                        checked,
                                        used_central,
                                        len(rel.cache),
                                    )

        if (central + 1) % 50 == 0:
            print(json.dumps({
                "progress_central": central + 1,
                "relation_rows_cached": len(rel.cache),
                "pattern_nodes_tested": checked,
            }, sort_keys=True), flush=True)

    return None, checked, used_central, len(rel.cache)


def pari_short_vectors(G, bound):
    pg = pari(G)
    try:
        result = pg.qfminim(ZZ(bound), None, 0)
    except TypeError:
        result = pari.qfminim(pg, ZZ(bound), None, 0)
    total = int(result[0])
    raw = Matrix(ZZ, result[2].sage())
    if raw.nrows() == G.nrows():
        pairs = raw.ncols()
    elif raw.ncols() == G.nrows():
        pairs = raw.nrows()
    else:
        raise RuntimeError("unexpected complement qfminim shape")
    return total, pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--central-limit", type=int, default=6660)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    G = SLOW.gram_matrix()
    result = {
        "status": "started",
        "algorithm": "PARI qfminim plus exact bitset Dynkin search",
        "candidate_lattice": "3.U4(3).2^2",
        "dimension": 30,
        "determinant": str(G.det()),
        "expected_determinant": str(3**15),
        "central_limit": args.central_limit,
        "truth_note": (
            "This is a lattice compatibility computation, not a Mordell-Weil "
            "realization or an elliptic-curve rank certificate."
        ),
    }
    if G.det() != 3**15:
        raise AssertionError("catalogue Gram determinant mismatch")

    vectors, total = enumerate_norm_six(G)
    result["norm6_vector_count"] = total
    result["norm6_vectors_up_to_sign"] = len(vectors)
    result["enumeration_seconds"] = time.time() - started

    mapping, checked, central_count, cache_size = search_embedding(
        vectors, G, args.central_limit
    )
    result.update({
        "pattern_nodes_tested": checked,
        "central_vectors_tried": central_count,
        "relation_rows_cached": cache_size,
        "embedding_found": mapping is not None,
    })

    if mapping is None:
        result["status"] = (
            "no_embedding" if central_count == len(vectors)
            else "no_embedding_in_tested_central_vectors"
        )
        result["elapsed_seconds"] = time.time() - started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({k: result.get(k) for k in (
            "status", "norm6_vectors_up_to_sign", "central_vectors_tried",
            "pattern_nodes_tested", "elapsed_seconds")}, sort_keys=True))
        return 0

    B = SLOW.orient_roots(mapping, vectors, G)
    target = SLOW.e8_cartan_scaled()
    subgram = B * G * B.transpose()
    if subgram != target:
        raise AssertionError("oriented sublattice is not E8(3)")

    primitive_index = B.index_in_saturation()
    K = SLOW.primitive_integer_kernel(B * G)
    if K.nrows() != 22 or K.ncols() != 30:
        raise AssertionError("orthogonal complement has wrong rank")
    complement_gram = K * G * K.transpose()
    block = B.stack(K)
    gluing_index = abs(block.det())
    determinant_identity_lhs = target.det() * complement_gram.det()
    determinant_identity_rhs = G.det() * gluing_index**2
    if determinant_identity_lhs != determinant_identity_rhs:
        raise AssertionError("determinant/gluing identity failed")

    complement_total, complement_pairs = pari_short_vectors(complement_gram, 6)
    result.update({
        "status": "pass",
        "embedding_indices_by_dynkin_node": {
            str(k): int(v) for k, v in mapping.items()
        },
        "embedding_basis_rows": [[int(x) for x in row] for row in B.rows()],
        "embedding_gram": [[int(x) for x in row] for row in subgram.rows()],
        "embedding_primitive_index": int(primitive_index),
        "embedding_is_primitive": bool(primitive_index == 1),
        "orthogonal_complement_rank": int(K.nrows()),
        "orthogonal_complement_basis_rows": [
            [int(x) for x in row] for row in K.rows()
        ],
        "orthogonal_complement_gram": [
            [int(x) for x in row] for row in complement_gram.rows()
        ],
        "orthogonal_complement_determinant": str(complement_gram.det()),
        "orthogonal_complement_smith_diagonal": [
            str(x) for x in complement_gram.smith_form()[0].diagonal()
        ],
        "orthogonal_complement_norm_at_most_6_vector_count": complement_total,
        "orthogonal_complement_norm_at_most_6_pairs_stored": complement_pairs,
        "direct_sum_gluing_index": str(gluing_index),
        "determinant_identity": {
            "det_E8_3_times_det_complement": str(determinant_identity_lhs),
            "det_L_times_gluing_index_squared": str(determinant_identity_rhs),
            "pass": True,
        },
        "elapsed_seconds": time.time() - started,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result.get(k) for k in (
        "status", "embedding_is_primitive", "orthogonal_complement_determinant",
        "direct_sum_gluing_index", "elapsed_seconds")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
