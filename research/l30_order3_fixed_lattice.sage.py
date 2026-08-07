#!/usr/bin/env sage-python
"""Construct an explicit order-three automorphism of the extremal
30-dimensional 3-modular lattice and certify its fixed lattice.

The cyclic-cubic rank-30 mechanism requires a deck transformation with
characteristic polynomial

    (x - 1)^8 (x^2 + x + 1)^11.

PARI ``qfauto`` computes exact lattice automorphisms from the catalogue Gram
matrix. Deterministic random words in these generators are projected to their
order-three parts. A successful element is checked entirely over ZZ, including
the fixed lattice, orthogonal complement, gluing index, and determinant
identities.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import random
import time
from pathlib import Path

from sage.all import Matrix, PolynomialRing, ZZ, identity_matrix, lcm, pari

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "l30_base", ROOT / "l30_e8_embedding.sage.py"
)
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)


def primitive_integer_kernel(A):
    Kq = A.right_kernel().basis_matrix()
    denominator = lcm([entry.denominator() for entry in Kq.list()] or [1])
    return (denominator * Kq).change_ring(ZZ).saturation()


def short_vector_data(G, bound=6):
    res = pari(G).qfminim(ZZ(bound), None, 0)
    total = int(res[0])
    raw = Matrix(ZZ, res[2].sage())
    pairs = raw.ncols() if raw.nrows() == G.nrows() else raw.nrows()
    return total, pairs


def random_word(rng, choices, identity, minimum=12, maximum=48):
    result = identity
    for _ in range(rng.randint(minimum, maximum)):
        result = result * choices[rng.randrange(len(choices))]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    G = BASE.gram_matrix()
    identity = identity_matrix(ZZ, 30)
    out = {
        "status": "started",
        "candidate_lattice": "3.U4(3).2^2",
        "dimension": 30,
        "determinant": str(G.det()),
        "expected_automorphism_group_order": "78382080",
        "requested_order_three_fixed_rank": 8,
        "trials": args.trials,
        "truth_note": (
            "An explicit lattice automorphism and fixed-lattice theorem do not "
            "by themselves realize this lattice as a Mordell-Weil lattice."
        ),
    }

    qfa = pari(G).qfauto()
    group_order = int(qfa[0])
    generators = [Matrix(ZZ, m.sage()) for m in qfa[1]]
    if group_order != 78382080:
        raise AssertionError(f"unexpected qfauto group order {group_order}")
    if not generators:
        raise AssertionError("qfauto returned no generators")
    for A in generators:
        if A.transpose() * G * A != G:
            raise AssertionError("qfauto generator does not preserve the Gram matrix")

    inverse_generators = []
    for A in generators:
        inverse = A.inverse()
        if any(entry.denominator() != 1 for entry in inverse.list()):
            raise AssertionError("a qfauto generator is not unimodular")
        inverse_generators.append(inverse.change_ring(ZZ))
    choices = generators + inverse_generators

    out["qfauto_seconds"] = time.time() - started
    out["automorphism_group_order"] = str(group_order)
    out["generator_count"] = len(generators)
    out["generator_determinants"] = [int(A.det()) for A in generators]

    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    wanted_charpoly = (x - 1)**8 * (x**2 + x + 1)**11
    chosen = None
    trial_records = []
    rng = random.Random(20260807)

    # If W lies in the automorphism group, then B=W^(|G|/3) satisfies B^3=1.
    # Nonidentity B therefore has exact order three. This avoids constructing
    # the full matrix group in GAP while preserving exactness.
    for trial in range(1, args.trials + 1):
        word = random_word(rng, choices, identity)
        candidate = word ** (group_order // 3)
        if candidate == identity:
            continue
        if candidate**3 != identity:
            raise AssertionError("the projected matrix does not have order three")
        if candidate.transpose() * G * candidate != G:
            raise AssertionError("the projected matrix is not a lattice isometry")
        fixed_rank = 30 - (candidate - identity).rank()
        cp = R(candidate.charpoly())
        trial_records.append({
            "trial": trial,
            "fixed_rank": int(fixed_rank),
            "characteristic_polynomial": str(cp),
        })
        if fixed_rank == 8 and cp == wanted_charpoly:
            chosen = candidate
            break

    out["order_three_elements_tested"] = len(trial_records)
    out["sampled_fixed_ranks"] = sorted({r["fixed_rank"] for r in trial_records})
    out["trial_records_tail"] = trial_records[-25:]

    if chosen is None:
        out["status"] = "no_matching_element_in_trials"
        out["elapsed_seconds"] = time.time() - started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
        print(json.dumps({k: out.get(k) for k in (
            "status", "automorphism_group_order", "generator_count",
            "order_three_elements_tested", "sampled_fixed_ranks",
            "elapsed_seconds")}, sort_keys=True))
        return 0

    fixed_basis = primitive_integer_kernel(chosen - identity)
    fixed_gram = fixed_basis * G * fixed_basis.transpose()
    fixed_target = BASE.e8_cartan_scaled()
    fixed_isom = pari(fixed_gram).qfisom(pari(fixed_target))
    if fixed_isom == 0:
        raise AssertionError("fixed lattice is not isometric to E8(3)")
    fixed_total, fixed_pairs = short_vector_data(fixed_gram, 6)

    complement_basis = primitive_integer_kernel(fixed_basis * G)
    complement_gram = complement_basis * G * complement_basis.transpose()
    complement_total, complement_pairs = short_vector_data(complement_gram, 6)
    block = fixed_basis.stack(complement_basis)
    gluing_index = abs(block.det())
    lhs = fixed_gram.det() * complement_gram.det()
    rhs = G.det() * gluing_index**2
    if lhs != rhs:
        raise AssertionError("fixed/complement gluing determinant identity failed")

    out.update({
        "status": "pass",
        "order_three_matrix_rows": [[int(v) for v in row] for row in chosen.rows()],
        "characteristic_polynomial": str(chosen.charpoly()),
        "fixed_rank": int(fixed_basis.nrows()),
        "fixed_basis_rows": [[int(v) for v in row] for row in fixed_basis.rows()],
        "fixed_gram": [[int(v) for v in row] for row in fixed_gram.rows()],
        "fixed_determinant": str(fixed_gram.det()),
        "fixed_isometric_to_E8_3": True,
        "fixed_isometry_matrix": [
            [int(v) for v in row] for row in Matrix(ZZ, fixed_isom.sage()).rows()
        ],
        "fixed_norm_at_most_6_vector_count": fixed_total,
        "fixed_norm_at_most_6_pairs_stored": fixed_pairs,
        "complement_rank": int(complement_basis.nrows()),
        "complement_basis_rows": [
            [int(v) for v in row] for row in complement_basis.rows()
        ],
        "complement_gram": [[int(v) for v in row] for row in complement_gram.rows()],
        "complement_determinant": str(complement_gram.det()),
        "complement_smith_diagonal": [
            str(v) for v in complement_gram.smith_form()[0].diagonal()
        ],
        "complement_norm_at_most_6_vector_count": complement_total,
        "complement_norm_at_most_6_pairs_stored": complement_pairs,
        "direct_sum_gluing_index": str(gluing_index),
        "determinant_identity": {
            "fixed_times_complement": str(lhs),
            "ambient_times_index_squared": str(rhs),
            "pass": True,
        },
        "elapsed_seconds": time.time() - started,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: out.get(k) for k in (
        "status", "fixed_rank", "fixed_determinant", "complement_determinant",
        "direct_sum_gluing_index", "elapsed_seconds")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
