#!/usr/bin/env sage-python
"""Exact lattice target for a cyclic-cubic rank-30 construction.

The lattice

    M0 = E8 (+) E8 (+) E8 (+) E6

has rank 30 and determinant 3.  An order-three isometry cyclically permutes
the three E8 factors and acts fixed-point-freely on E6 via the fourth power
of an E6 Coxeter element.  The fixed lattice is E8(3), exactly the lattice
forced by pulling a split E8 rational elliptic surface through a cubic cover.

The script also gives an explicit discriminant-form gluing showing that

    NS0 = <O,F> (+) M0(-1),  O^2=-3, O.F=1, F^2=0,

with transcendental candidate E6, embeds primitively in an odd unimodular
lattice of signature (7,31), the topological H^2 lattice required for a
chi=3 elliptic surface over a genus-one base.

This is an inverse-lattice feasibility theorem, not an elliptic-surface or
rank-30 realization theorem.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import (
    CartanMatrix,
    Matrix,
    PolynomialRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    pari,
    vector,
)


def primitive_integer_kernel(A):
    Kq = A.right_kernel().basis_matrix()
    denominator = lcm([entry.denominator() for entry in Kq.list()] or [1])
    return (denominator * Kq).change_ring(ZZ).saturation()


def simple_reflection(C, i):
    """Simple-root reflection acting on column coordinates."""
    n = C.nrows()
    S = identity_matrix(ZZ, n)
    for j in range(n):
        S[i, j] -= C[i, j]
    if S.transpose() * C * S != C:
        raise AssertionError("simple reflection does not preserve the Cartan form")
    return S


def coxeter_matrix(C):
    out = identity_matrix(ZZ, C.nrows())
    for i in range(C.nrows()):
        out = out * simple_reflection(C, i)
    return out


def matrix_signature(G):
    values = G.change_ring(QQ).eigenvalues()
    positive = sum(1 for x in values if x > 0)
    negative = sum(1 for x in values if x < 0)
    zero = len(values) - positive - negative
    return [positive, negative, zero]


def qfminim_count(G, bound):
    result = pari(G).qfminim(ZZ(bound), None, 0)
    total = int(result[0])
    raw = Matrix(ZZ, result[2].sage())
    if raw.nrows() == G.nrows():
        pairs = raw.ncols()
    elif raw.ncols() == G.nrows():
        pairs = raw.nrows()
    else:
        raise RuntimeError(f"unexpected qfminim shape {raw.dimensions()}")
    return total, pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    E8 = Matrix(ZZ, CartanMatrix(["E", 8]))
    E6 = Matrix(ZZ, CartanMatrix(["E", 6]))
    M = block_diagonal_matrix(E8, E8, E8, E6)
    if M.nrows() != 30 or M.det() != 3:
        raise AssertionError("target lattice invariants are wrong")

    I8 = identity_matrix(ZZ, 8)
    Z8 = Matrix(ZZ, 8, 8, 0)
    Z86 = Matrix(ZZ, 8, 6, 0)
    Z68 = Matrix(ZZ, 6, 8, 0)

    coxeter = coxeter_matrix(E6)
    sigma6 = coxeter**4
    if sigma6**3 != identity_matrix(ZZ, 6) or sigma6 == identity_matrix(ZZ, 6):
        raise AssertionError("the E6 Coxeter fourth power is not order three")
    if sigma6.transpose() * E6 * sigma6 != E6:
        raise AssertionError("the E6 order-three matrix is not an isometry")
    if (sigma6 - identity_matrix(ZZ, 6)).rank() != 6:
        raise AssertionError("the E6 order-three isometry has fixed vectors")

    sigma = Matrix(ZZ, 30, 30, 0)
    # Column action (v1,v2,v3,w) -> (v3,v1,v2,sigma6*w).
    sigma[0:8, 16:24] = I8
    sigma[8:16, 0:8] = I8
    sigma[16:24, 8:16] = I8
    sigma[24:30, 24:30] = sigma6

    I30 = identity_matrix(ZZ, 30)
    if sigma**3 != I30 or sigma == I30:
        raise AssertionError("the target automorphism is not order three")
    if sigma.transpose() * M * sigma != M:
        raise AssertionError("the target automorphism does not preserve M0")

    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    expected_cp = (x - 1)**8 * (x**2 + x + 1)**11
    if R(sigma.charpoly()) != expected_cp:
        raise AssertionError("the cyclic-cubic characteristic polynomial is wrong")

    fixed_basis = primitive_integer_kernel(sigma - I30)
    fixed_gram = fixed_basis * M * fixed_basis.transpose()
    fixed_target = 3 * E8
    fixed_isom = pari(fixed_gram).qfisom(pari(fixed_target))
    if fixed_isom == 0:
        raise AssertionError("the fixed lattice is not E8(3)")

    complement_basis = primitive_integer_kernel(fixed_basis * M)
    complement_gram = complement_basis * M * complement_basis.transpose()
    block = fixed_basis.stack(complement_basis)
    direct_sum_index = abs(block.det())
    if fixed_gram.det() != 3**8:
        raise AssertionError("fixed-lattice determinant is not 3^8")
    if complement_gram.det() != 3**9:
        raise AssertionError("trace-zero determinant is not 3^9")
    if direct_sum_index != 3**8:
        raise AssertionError("fixed/complement gluing index is not 3^8")
    if fixed_gram.det() * complement_gram.det() != M.det() * direct_sum_index**2:
        raise AssertionError("fixed/complement determinant identity failed")

    ambient_roots, ambient_pairs = qfminim_count(M, 2)
    fixed_minimal, fixed_pairs = qfminim_count(fixed_gram, 6)
    complement_roots, complement_pairs = qfminim_count(complement_gram, 2)
    if ambient_roots != 792 or ambient_pairs != 396:
        raise AssertionError("M0 must have exactly 792 roots")
    if fixed_minimal != 240 or fixed_pairs != 120:
        raise AssertionError("E8(3) must have exactly 240 minimal vectors")
    if complement_roots != 72 or complement_pairs != 36:
        raise AssertionError("the trace-zero lattice must have the 72 E6 roots")

    # The section-fibre lattice for chi=3 is odd unimodular of signature (1,1).
    OF = Matrix(ZZ, [[-3, 1], [1, 0]])
    NS = block_diagonal_matrix(OF, -M)
    T = E6
    H0 = block_diagonal_matrix(NS, T)
    if NS.det() != -3 or T.det() != 3 or H0.det() != -9:
        raise AssertionError("NS/transcendental discriminants are wrong")

    # Glue a generator of A_{E6(-1)} to the same fundamental-weight class in
    # A_{E6(+1)}.  The two discriminant quadratic values cancel.
    e1 = vector(ZZ, [1, 0, 0, 0, 0, 0])
    weight = E6.inverse() * e1
    if any((3 * a).denominator() != 1 for a in weight):
        raise AssertionError("the selected E6 weight does not have order three")
    if all(a.denominator() == 1 for a in weight):
        raise AssertionError("the selected E6 weight is integral")

    glue = vector(QQ, [0] * 38)
    negative_e6_start = 2 + 24
    positive_e6_start = 32
    for i in range(6):
        glue[negative_e6_start + i] = weight[i]
        glue[positive_e6_start + i] = weight[i]
    if any(a.denominator() != 1 for a in 3 * glue):
        raise AssertionError("the gluing vector does not have order three")
    pairings = glue * H0
    if any(a.denominator() != 1 for a in pairings):
        raise AssertionError("the gluing vector is not in the dual lattice")
    glue_norm = (glue * H0).dot_product(glue)
    if glue_norm != 0:
        raise AssertionError("the discriminant-form gluing vector is not isotropic")

    nonintegral_coordinate = next(
        i for i, value in enumerate(glue) if value.denominator() == 3
    )
    super_basis_rows = [glue]
    for i in range(38):
        if i != nonintegral_coordinate:
            row = vector(QQ, [0] * 38)
            row[i] = 1
            super_basis_rows.append(row)
    super_basis = Matrix(QQ, super_basis_rows)
    if abs(super_basis.det()) != QQ(1, 3):
        raise AssertionError("the proposed overlattice does not have index three")
    H = super_basis * H0 * super_basis.transpose()
    if any(value.denominator() != 1 for value in H.list()):
        raise AssertionError("the glued intersection form is not integral")
    H = H.change_ring(ZZ)
    if abs(H.det()) != 1:
        raise AssertionError("the glued intersection form is not unimodular")
    signature = matrix_signature(H)
    if signature != [7, 31, 0]:
        raise AssertionError(f"unexpected H2 signature {signature}")
    if all(H[i, i] % 2 == 0 for i in range(38)):
        raise AssertionError("the glued H2 lattice should be odd")

    result = {
        "status": "pass",
        "theorem": "explicit cyclic-cubic rank-30 target lattice",
        "mordell_weil_candidate": "E8^3 direct-sum E6",
        "rank": 30,
        "determinant": "3",
        "minimum": 2,
        "root_count": ambient_roots,
        "order_three_matrix_rows": [[int(v) for v in row] for row in sigma.rows()],
        "order_three_characteristic_polynomial": str(sigma.charpoly()),
        "fixed_rank": fixed_basis.nrows(),
        "fixed_basis_rows": [[int(v) for v in row] for row in fixed_basis.rows()],
        "fixed_gram": [[int(v) for v in row] for row in fixed_gram.rows()],
        "fixed_isometric_to_E8_3": True,
        "fixed_determinant": str(fixed_gram.det()),
        "fixed_minimal_vector_count": fixed_minimal,
        "trace_zero_rank": complement_basis.nrows(),
        "trace_zero_basis_rows": [
            [int(v) for v in row] for row in complement_basis.rows()
        ],
        "trace_zero_gram": [
            [int(v) for v in row] for row in complement_gram.rows()
        ],
        "trace_zero_determinant": str(complement_gram.det()),
        "trace_zero_root_count": complement_roots,
        "fixed_plus_trace_zero_index": str(direct_sum_index),
        "section_orbits_of_roots": 240 + 24,
        "neron_severi_candidate": {
            "rank": NS.nrows(),
            "signature": matrix_signature(NS),
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
            "glue_vector": [str(v) for v in glue],
            "glue_norm": str(glue_norm),
            "overlattice_index": "3",
            "unimodular_determinant": str(H.det()),
            "signature": signature,
            "odd": True,
        },
        "truth_note": (
            "All lattice and topological compatibility checks are exact. The "
            "remaining problem is algebraic-geometric realization as the "
            "Mordell-Weil lattice of a cyclic cubic pullback, followed by rational "
            "specialization and a 30-point independence certificate."
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
        "H2_signature": signature,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
