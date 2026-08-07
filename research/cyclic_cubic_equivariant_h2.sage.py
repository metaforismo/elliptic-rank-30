#!/usr/bin/env sage-python
"""Exact equivariant H^2 compatibility for the cyclic-cubic rank-30 target.

Starting from the target Mordell--Weil lattice

    M0 = E8^3 (+) E6,

this script constructs:

* the order-three action with fixed Mordell--Weil lattice E8(3);
* the candidate Neron--Severi lattice
      [[-3,1],[1,0]] (+) M0(-1);
* the transcendental candidate E6;
* the index-three discriminant gluing to an odd unimodular lattice of
  signature (7,31);
* an integral order-three isometry of the glued H^2 lattice.

The resulting H^2 characteristic polynomial is

    (x-1)^10 (x^2+x+1)^14,

so Tr(sigma|H^2)=-4.  Together with traces -1 on H^1 and H^3, the
topological Lefschetz number is zero, exactly matching three fixed smooth
elliptic branch fibres.

This is an equivariant lattice/topology theorem.  It does not construct the
algebraic surface or 30 rational points.
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
    vector,
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


def primitive_integer_kernel(A):
    Kq = A.right_kernel().basis_matrix()
    denominator = lcm([entry.denominator() for entry in Kq.list()] or [1])
    return (denominator * Kq).change_ring(ZZ).saturation()


def exact_signature(G):
    eigenvalues = G.change_ring(QQ).eigenvalues()
    positive = sum(value > 0 for value in eigenvalues)
    negative = sum(value < 0 for value in eigenvalues)
    return [positive, negative, len(eigenvalues) - positive - negative]


def target_data():
    E8 = Matrix(ZZ, CartanMatrix(["E", 8]))
    E6 = Matrix(ZZ, CartanMatrix(["E", 6]))
    M = block_diagonal_matrix(E8, E8, E8, E6)

    sigma6 = coxeter(E6)**4
    I6 = identity_matrix(ZZ, 6)
    assert sigma6**3 == I6 and sigma6 != I6
    assert sigma6.transpose() * E6 * sigma6 == E6
    assert (sigma6 - I6).rank() == 6

    I8 = identity_matrix(ZZ, 8)
    sigma_mw = Matrix(ZZ, 30, 30, 0)
    sigma_mw[0:8, 16:24] = I8
    sigma_mw[8:16, 0:8] = I8
    sigma_mw[16:24, 8:16] = I8
    sigma_mw[24:30, 24:30] = sigma6
    I30 = identity_matrix(ZZ, 30)
    assert sigma_mw**3 == I30 and sigma_mw != I30
    assert sigma_mw.transpose() * M * sigma_mw == M

    OF = Matrix(ZZ, [[-3, 1], [1, 0]])
    NS = block_diagonal_matrix(OF, -M)
    T = E6
    H0 = block_diagonal_matrix(NS, T)

    sigma_h0 = block_diagonal_matrix(
        identity_matrix(ZZ, 2), sigma_mw, sigma6
    )
    I38 = identity_matrix(ZZ, 38)
    assert sigma_h0**3 == I38 and sigma_h0 != I38
    assert sigma_h0.transpose() * H0 * sigma_h0 == H0

    # Equal order-three fundamental-weight classes in E6(-1) and E6(+1).
    weight = E6.inverse() * vector(ZZ, [1, 0, 0, 0, 0, 0])
    assert all((3 * value).denominator() == 1 for value in weight)
    assert any(value.denominator() == 3 for value in weight)

    glue_class = vector(QQ, [0] * 38)
    for i in range(6):
        glue_class[26 + i] = weight[i]
        glue_class[32 + i] = weight[i]
    assert all(value.denominator() == 1 for value in glue_class * H0)
    assert (glue_class * H0).dot_product(glue_class) == 0

    # Normalize one coefficient to +/-1/3 by adding an integral basis vector.
    pivot = next(i for i, value in enumerate(glue_class) if value.denominator() == 3)
    glue_basis = vector(QQ, glue_class)
    numerator = glue_basis[pivot].numerator()
    desired = 1 if numerator % 3 == 1 else -1
    glue_basis[pivot] += QQ(desired - numerator) / 3
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
    assert H.det() == -1
    assert exact_signature(H) == [7, 31, 0]
    assert any(H[i, i] % 2 for i in range(38))

    # Row-coordinate action in the glued basis.
    action_q = change * sigma_h0.transpose() * change.inverse()
    assert all(value.denominator() == 1 for value in action_q.list())
    action = action_q.change_ring(ZZ)
    assert action**3 == I38 and action != I38
    assert action * H * action.transpose() == H

    return {
        "E8": E8,
        "E6": E6,
        "M": M,
        "sigma6": sigma6,
        "sigma_mw": sigma_mw,
        "OF": OF,
        "NS": NS,
        "T": T,
        "H0": H0,
        "sigma_h0": sigma_h0,
        "glue_class": glue_class,
        "glue_basis": glue_basis,
        "change": change,
        "H": H,
        "action": action,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = target_data()
    H = data["H"]
    action = data["action"]
    sigma_mw = data["sigma_mw"]

    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    mw_expected = (x - 1)**8 * (x**2 + x + 1)**11
    h2_expected = (x - 1)**10 * (x**2 + x + 1)**14
    assert R(sigma_mw.charpoly()) == mw_expected
    assert R(action.charpoly()) == h2_expected

    fixed_h2 = primitive_integer_kernel(action - identity_matrix(ZZ, 38))
    fixed_h2_gram = fixed_h2 * H * fixed_h2.transpose()
    assert fixed_h2.nrows() == 10
    assert exact_signature(fixed_h2_gram) == [1, 9, 0]

    trace_h2 = int(action.trace())
    assert trace_h2 == -4
    trace_h0 = 1
    trace_h1 = -1  # quotient of the genus-one base by C3 has genus zero
    trace_h3 = -1
    trace_h4 = 1
    lefschetz = trace_h0 - trace_h1 + trace_h2 - trace_h3 + trace_h4
    assert lefschetz == 0

    result = {
        "status": "pass",
        "theorem": "equivariant H2 compatibility for the cyclic-cubic rank-30 target",
        "mordell_weil_characteristic_polynomial_factored": "(x-1)^8*(x^2+x+1)^11",
        "H2_characteristic_polynomial_factored": "(x-1)^10*(x^2+x+1)^14",
        "H2_rank": 38,
        "H2_signature": exact_signature(H),
        "H2_determinant": str(H.det()),
        "H2_odd_unimodular": True,
        "H2_order_three_matrix_rows": [[int(v) for v in row] for row in action.rows()],
        "H2_fixed_rank": 10,
        "H2_fixed_signature": exact_signature(fixed_h2_gram),
        "H2_fixed_gram": [[int(v) for v in row] for row in fixed_h2_gram.rows()],
        "H2_trace": trace_h2,
        "cohomology_traces": {
            "H0": trace_h0,
            "H1": trace_h1,
            "H2": trace_h2,
            "H3": trace_h3,
            "H4": trace_h4,
        },
        "topological_lefschetz_number": lefschetz,
        "expected_fixed_locus": {
            "description": "three smooth elliptic branch fibres",
            "euler_characteristic": 0,
            "matches_lefschetz_number": True,
        },
        "glue_class_stable": True,
        "truth_note": (
            "This proves exact equivariant lattice and topological compatibility. "
            "It does not prove that an algebraic elliptic surface realizes the "
            "data, nor does it exhibit 30 independent rational points."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "H2_fixed_rank": result["H2_fixed_rank"],
        "H2_trace": result["H2_trace"],
        "topological_lefschetz_number": result["topological_lefschetz_number"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
