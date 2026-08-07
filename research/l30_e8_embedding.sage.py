#!/usr/bin/env sage-python
"""Exact lattice compatibility test for the cyclic-cubic rank-30 route.

The candidate lattice is the 30-dimensional integral 3-modular lattice
catalogued as 3.U4(3).2^2.  A cyclic cubic pullback of a split E8 rational
elliptic surface necessarily contains the pulled-back generic lattice E8(3).
This script asks the decisive finite question: does the candidate lattice
contain a primitive sublattice with Gram matrix 3*Cartan(E8)?

The computation is certificate-producing.  It verifies the catalogue
invariants, enumerates every norm-6 vector, performs an exact Dynkin-pattern
search, orients the roots, checks primitivity, and computes the orthogonal
complement and gluing index.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from sage.all import IntegralLattice, Matrix, ZZ, identity_matrix, lcm, vector

LOWER = [
    [6],
    [3, 6],
    [2, 3, 6],
    [2, 3, 3, 6],
    [2, 3, 2, 2, 6],
    [3, 3, 2, 2, 2, 6],
    [2, 3, 2, 2, 2, 2, 6],
    [2, 3, 0, 3, 2, 2, 2, 6],
    [2, 3, 2, 2, 2, 2, 3, 2, 6],
    [2, 1, 2, 0, 2, 0, 1, 0, 2, 6],
    [0, 1, 2, 2, 2, 0, 1, 0, 1, 2, 6],
    [2, 3, 2, 2, 2, 2, 3, 2, 3, 0, 1, 6],
    [2, 3, 2, 2, 3, 2, 2, 2, 2, 2, 1, 2, 6],
    [3, 3, 2, 2, 2, 0, 2, 2, 2, 2, 1, 2, 2, 6],
    [2, 3, 2, 2, 2, 2, 3, 2, 0, 1, 2, 0, 2, 2, 6],
    [2, 2, 2, 3, 2, 1, 2, 1, 1, 2, 2, 0, 0, 2, 2, 6],
    [2, 3, 3, 3, 2, 2, 2, 0, 2, 1, 2, 2, 2, 2, 2, 1, 6],
    [0, 2, 2, 1, 2, 0, 2, 1, 1, 1, 1, 0, 0, 2, 2, 3, 0, 6],
    [1, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 2, 0, 0, 1, 2, 2, 2, 6],
    [2, 1, 3, 2, 2, 2, 1, 1, 1, 2, 2, 1, 2, 0, 2, 2, 1, 1, 1, 6],
    [0, 2, 2, 1, 2, 2, 2, 1, 2, 0, 2, 2, 0, 0, 1, 2, 0, 2, 2, 2, 6],
    [0, 1, 1, 2, 0, 0, 3, 2, 2, 1, 1, 1, 0, 1, 2, 2, 0, 3, 2, 1, 1, 6],
    [2, 2, 2, 3, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 3, 1, 2, 2, 3, 2, 2, 6],
    [0, 1, 1, 1, 2, 2, 2, 2, 0, -1, 2, 1, 1, -1, 2, 0, 1, 2, 2, 2, 2, 2, 0, 6],
    [2, 3, 2, 2, 0, 2, 2, 2, 2, 1, 1, 2, 3, 2, 2, 0, 2, 0, 0, 1, 0, 2, 2, 1, 6],
    [3, 3, 2, 2, 2, 3, 2, 2, 2, 1, 0, 2, 2, 0, 2, 0, 2, 0, 1, 2, 1, 0, 2, 1, 2, 6],
    [1, 1, 3, 2, 2, 2, 1, 1, 2, 2, 2, 0, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 0, 6],
    [2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 1, 0, 2, 2, 2, 2, 2, 1, 0, 2, 2, 0, 2, 2, 6],
    [3, 2, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 0, 2, 2, 3, 2, 2, 3, 2, 2, 2, 1, 2, 6],
    [2, 2, 2, 1, 0, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 0, 1, 2, 2, 2, 2, 2, 0, 1, 0, 1, 1, 3, 6],
]

# E8 Dynkin diagram: 0-2-3-4-5-6-7 with 1-2.
ADJ = {
    0: {2},
    1: {2},
    2: {0, 1, 3},
    3: {2, 4},
    4: {3, 5},
    5: {4, 6},
    6: {5, 7},
    7: {6},
}
SEARCH_ORDER = [2, 0, 1, 3, 4, 5, 6, 7]


def gram_matrix():
    G = Matrix(ZZ, 30, 30)
    for i, row in enumerate(LOWER):
        for j, value in enumerate(row):
            G[i, j] = value
            G[j, i] = value
    return G


def e8_cartan_scaled():
    C = 6 * identity_matrix(ZZ, 8)
    for i in range(8):
        for j in ADJ[i]:
            C[i, j] = -3
    return C


def search_embedding(short_vectors, G, central_limit):
    gv = [v * G for v in short_vectors]

    def ip(i, j):
        return gv[i].dot_product(short_vectors[j])

    checked = 0
    for central_index in range(min(central_limit, len(short_vectors))):
        assigned = {2: central_index}
        used = {central_index}

        def compatible(node, candidate):
            for other_node, other_candidate in assigned.items():
                value = ip(candidate, other_candidate)
                if other_node in ADJ[node]:
                    if abs(value) != 3:
                        return False
                elif value != 0:
                    return False
            return True

        def recurse(position):
            nonlocal checked
            if position == len(SEARCH_ORDER):
                return dict(assigned)
            node = SEARCH_ORDER[position]
            for candidate in range(len(short_vectors)):
                if candidate in used:
                    continue
                checked += 1
                if not compatible(node, candidate):
                    continue
                # Remove the symmetry interchanging the two short arms.
                if node == 1 and assigned.get(0, -1) >= candidate:
                    continue
                assigned[node] = candidate
                used.add(candidate)
                answer = recurse(position + 1)
                if answer is not None:
                    return answer
                used.remove(candidate)
                del assigned[node]
            return None

        answer = recurse(1)
        if answer is not None:
            return answer, checked, central_index + 1
    return None, checked, min(central_limit, len(short_vectors))


def orient_roots(mapping, short_vectors, G):
    # Root 2 is the starting vertex. Orient along the tree so every edge is -3.
    oriented = {2: short_vectors[mapping[2]]}
    queue = [2]
    while queue:
        parent = queue.pop(0)
        for child in sorted(ADJ[parent]):
            if child in oriented:
                continue
            raw = short_vectors[mapping[child]]
            value = (oriented[parent] * G).dot_product(raw)
            if abs(value) != 3:
                raise AssertionError("found pattern lost an E8 edge")
            oriented[child] = raw if value == -3 else -raw
            queue.append(child)
    return Matrix(ZZ, [oriented[i] for i in range(8)])


def primitive_integer_kernel(A):
    Kq = A.right_kernel().basis_matrix()
    denominator = lcm([entry.denominator() for entry in Kq.list()] or [1])
    Ki = (denominator * Kq).change_ring(ZZ)
    return Ki.saturation()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--central-limit", type=int, default=64)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    G = gram_matrix()
    result = {
        "status": "started",
        "candidate_lattice": "3.U4(3).2^2",
        "source": "Nebe-Sloane Catalogue of Lattices, neb26",
        "dimension": 30,
        "determinant": str(G.det()),
        "expected_determinant": str(3**15),
        "smith_diagonal": [str(x) for x in G.smith_form()[0].diagonal()],
        "central_limit": args.central_limit,
        "truth_note": (
            "This is a lattice-compatibility computation.  It is not an elliptic "
            "curve, a Mordell-Weil realization theorem, or a rank-30 certificate."
        ),
    }
    if G.det() != 3**15:
        raise AssertionError("catalogue Gram determinant mismatch")

    L = IntegralLattice(G)
    short_by_norm = L.short_vectors(7, up_to_sign_flag=True)
    short_vectors = [vector(ZZ, v) for v in short_by_norm[6]]
    result["norm6_vectors_up_to_sign"] = len(short_vectors)
    result["kissing_number"] = 2 * len(short_vectors)
    if len(short_vectors) != 6660:
        raise AssertionError("catalogue kissing number mismatch")

    mapping, checked, central_count = search_embedding(short_vectors, G, args.central_limit)
    result["pattern_candidates_checked"] = checked
    result["central_vectors_tried"] = central_count
    result["embedding_found"] = mapping is not None

    if mapping is None:
        result["status"] = "no_embedding_in_tested_central_orbits"
        result["elapsed_seconds"] = time.time() - started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({k: result.get(k) for k in (
            "status", "norm6_vectors_up_to_sign", "central_vectors_tried",
            "pattern_candidates_checked")}, sort_keys=True))
        return 0

    B = orient_roots(mapping, short_vectors, G)
    target = e8_cartan_scaled()
    subgram = B * G * B.transpose()
    if subgram != target:
        raise AssertionError("oriented sublattice is not E8(3)")

    primitive_index = B.index_in_saturation()
    K = primitive_integer_kernel(B * G)
    if K.nrows() != 22 or K.ncols() != 30:
        raise AssertionError("orthogonal complement has wrong rank")
    complement_gram = K * G * K.transpose()
    block = B.stack(K)
    gluing_index = abs(block.det())
    determinant_identity_lhs = target.det() * complement_gram.det()
    determinant_identity_rhs = G.det() * gluing_index**2
    if determinant_identity_lhs != determinant_identity_rhs:
        raise AssertionError("determinant/gluing identity failed")

    complement_lattice = IntegralLattice(complement_gram)
    complement_short = complement_lattice.short_vectors(7, up_to_sign_flag=True)

    result.update({
        "status": "pass",
        "embedding_indices_by_dynkin_node": {str(k): int(v) for k, v in mapping.items()},
        "embedding_basis_rows": [[int(x) for x in row] for row in B.rows()],
        "embedding_gram": [[int(x) for x in row] for row in subgram.rows()],
        "embedding_primitive_index": int(primitive_index),
        "embedding_is_primitive": bool(primitive_index == 1),
        "orthogonal_complement_rank": int(K.nrows()),
        "orthogonal_complement_basis_rows": [[int(x) for x in row] for row in K.rows()],
        "orthogonal_complement_gram": [[int(x) for x in row] for row in complement_gram.rows()],
        "orthogonal_complement_determinant": str(complement_gram.det()),
        "orthogonal_complement_smith_diagonal": [str(x) for x in complement_gram.smith_form()[0].diagonal()],
        "orthogonal_complement_norm6_vectors_up_to_sign": len(complement_short[6]),
        "orthogonal_complement_kissing_number_at_6": 2 * len(complement_short[6]),
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
