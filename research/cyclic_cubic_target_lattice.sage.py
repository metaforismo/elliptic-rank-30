#!/usr/bin/env sage-python
"""Exact inverse-lattice target for a cyclic-cubic rank-30 construction.

We certify the rank-30 lattice M0 = E8^3 (+) E6, an order-three isometry
with fixed lattice E8(3), its rank-22 trace-zero complement, and an explicit
discriminant-form gluing into the odd unimodular H^2 lattice of signature
(7,31) required for an elliptic surface with chi=3 over a genus-one base.

This is a lattice/topology feasibility theorem, not an elliptic-surface or
rank-30 realization theorem.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import (
    CartanMatrix, Matrix, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, pari, vector,
)


def reflection(C, i):
    S = identity_matrix(ZZ, C.nrows())
    for j in range(C.ncols()):
        S[i, j] -= C[i, j]
    assert S.transpose() * C * S == C
    return S


def coxeter(C):
    W = identity_matrix(ZZ, C.nrows())
    for i in range(C.nrows()):
        W *= reflection(C, i)
    return W


def signature(G):
    eigenvalues = G.change_ring(QQ).eigenvalues()
    positive = sum(value > 0 for value in eigenvalues)
    negative = sum(value < 0 for value in eigenvalues)
    return [positive, negative, len(eigenvalues) - positive - negative]


def qf_count(G, bound):
    result = pari(G).qfminim(ZZ(bound), None, 0)
    total = int(result[0])
    raw = Matrix(ZZ, result[2].sage())
    pairs = raw.ncols() if raw.nrows() == G.nrows() else raw.nrows()
    return total, pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    E8 = Matrix(ZZ, CartanMatrix(["E", 8]))
    E6 = Matrix(ZZ, CartanMatrix(["E", 6]))
    M = block_diagonal_matrix(E8, E8, E8, E6)
    assert M.nrows() == 30 and M.det() == 3

    # Fixed-point-free order-three Weyl isometry on E6.
    sigma6 = coxeter(E6)**4
    I6 = identity_matrix(ZZ, 6)
    assert sigma6**3 == I6 and sigma6 != I6
    assert sigma6.transpose() * E6 * sigma6 == E6
    assert (sigma6 - I6).rank() == 6

    # Cycle the three E8 factors.
    I8 = identity_matrix(ZZ, 8)
    sigma = Matrix(ZZ, 30, 30, 0)
    sigma[0:8, 16:24] = I8
    sigma[8:16, 0:8] = I8
    sigma[16:24, 8:16] = I8
    sigma[24:30, 24:30] = sigma6
    I30 = identity_matrix(ZZ, 30)
    assert sigma**3 == I30 and sigma != I30
    assert sigma.transpose() * M * sigma == M

    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    expected_cp = (x - 1)**8 * (x**2 + x + 1)**11
    assert R(sigma.charpoly()) == expected_cp

    # Fixed lattice (v,v,v,0) = E8(3).
    fixed = Matrix(ZZ, 8, 30, 0)
    fixed[:, 0:8] = I8
    fixed[:, 8:16] = I8
    fixed[:, 16:24] = I8
    assert fixed * sigma.transpose() == fixed
    assert fixed.index_in_saturation() == 1
    fixed_gram = fixed * M * fixed.transpose()
    assert fixed_gram == 3 * E8 and fixed_gram.det() == 3**8

    # Orthogonal trace-zero lattice: E8 tensor A2, plus E6.
    trace = Matrix(ZZ, 22, 30, 0)
    trace[0:8, 0:8] = I8
    trace[0:8, 8:16] = -I8
    trace[8:16, 8:16] = I8
    trace[8:16, 16:24] = -I8
    trace[16:22, 24:30] = I6
    assert fixed * M * trace.transpose() == 0
    assert trace.index_in_saturation() == 1
    trace_gram = trace * M * trace.transpose()
    assert trace_gram.det() == 3**9

    glue_index = abs(fixed.stack(trace).det())
    assert glue_index == 3**8
    assert fixed_gram.det() * trace_gram.det() == M.det() * glue_index**2

    ambient_roots, ambient_pairs = qf_count(M, 2)
    fixed_minimal, fixed_pairs = qf_count(fixed_gram, 6)
    trace_roots, trace_pairs = qf_count(trace_gram, 2)
    assert (ambient_roots, ambient_pairs) == (792, 396)
    assert (fixed_minimal, fixed_pairs) == (240, 120)
    assert (trace_roots, trace_pairs) == (72, 36)

    # Candidate Neron-Severi lattice and transcendental complement.
    OF = Matrix(ZZ, [[-3, 1], [1, 0]])
    NS = block_diagonal_matrix(OF, -M)
    T = E6
    H0 = block_diagonal_matrix(NS, T)
    assert NS.det() == -3 and T.det() == 3 and H0.det() == -9

    # Glue equal fundamental-weight classes in E6(-1) and E6(+1).
    weight = E6.inverse() * vector(ZZ, [1, 0, 0, 0, 0, 0])
    assert all((3 * value).denominator() == 1 for value in weight)
    assert any(value.denominator() == 3 for value in weight)

    glue = vector(QQ, [0] * 38)
    for i in range(6):
        glue[26 + i] = weight[i]
        glue[32 + i] = weight[i]
    assert all(value.denominator() == 1 for value in 3 * glue)
    assert all(value.denominator() == 1 for value in glue * H0)
    assert (glue * H0).dot_product(glue) == 0

    # Choose an equivalent representative whose pivot coefficient is +/-1/3.
    pivot = next(i for i, value in enumerate(glue) if value.denominator() == 3)
    glue_basis = vector(QQ, glue)
    numerator = glue_basis[pivot].numerator()
    target_numerator = 1 if numerator % 3 == 1 else -1
    glue_basis[pivot] += QQ(target_numerator - numerator) / 3
    assert abs(glue_basis[pivot]) == QQ(1) / 3
    assert all(value.denominator() == 1 for value in glue_basis * H0)
    assert (glue_basis * H0).dot_product(glue_basis).denominator() == 1

    rows = [glue_basis]
    for i in range(38):
        if i == pivot:
            continue
        row = vector(QQ, [0] * 38)
        row[i] = 1
        rows.append(row)
    change = Matrix(QQ, rows)
    assert abs(change.det()) == QQ(1) / 3
    H = change * H0 * change.transpose()
    assert all(value.denominator() == 1 for value in H.list())
    H = H.change_ring(ZZ)
    assert abs(H.det()) == 1
    H_signature = signature(H)
    assert H_signature == [7, 31, 0]
    assert any(H[i, i] % 2 for i in range(38))

    result = {
        "status": "pass",
        "theorem": "explicit cyclic-cubic rank-30 target lattice",
        "mordell_weil_candidate": "E8^3 direct-sum E6",
        "rank": 30,
        "determinant": "3",
        "minimum": 2,
        "root_count": ambient_roots,
        "order_three_characteristic_polynomial": str(sigma.charpoly()),
        "order_three_matrix_rows": [[int(v) for v in row] for row in sigma.rows()],
        "fixed_rank": 8,
        "fixed_gram": [[int(v) for v in row] for row in fixed_gram.rows()],
        "fixed_isometric_to_E8_3": True,
        "fixed_determinant": str(fixed_gram.det()),
        "fixed_minimal_vector_count": fixed_minimal,
        "trace_zero_rank": 22,
        "trace_zero_gram": [[int(v) for v in row] for row in trace_gram.rows()],
        "trace_zero_determinant": str(trace_gram.det()),
        "trace_zero_root_count": trace_roots,
        "fixed_plus_trace_zero_index": str(glue_index),
        "section_orbits_of_roots": 264,
        "neron_severi_candidate": {
            "rank": 32,
            "signature": signature(NS),
            "determinant": str(NS.det()),
            "trivial_lattice_gram": [[-3, 1], [1, 0]],
        },
        "transcendental_candidate": {
            "lattice": "E6",
            "rank": 6,
            "determinant": str(T.det()),
        },
        "topological_gluing": {
            "unglued_determinant": str(H0.det()),
            "discriminant_glue_class": [str(v) for v in glue],
            "basis_representative": [str(v) for v in glue_basis],
            "overlattice_index": "3",
            "unimodular_determinant": str(H.det()),
            "signature": H_signature,
            "odd": True,
        },
        "truth_note": (
            "All lattice and topological checks are exact. Remaining: realize "
            "this data by a cyclic cubic elliptic surface over Q, construct its "
            "30 sections, specialize, and certify independence."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "rank": result["rank"],
        "fixed_rank": result["fixed_rank"],
        "trace_zero_rank": result["trace_zero_rank"],
        "fixed_plus_trace_zero_index": result["fixed_plus_trace_zero_index"],
        "H2_signature": H_signature,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
