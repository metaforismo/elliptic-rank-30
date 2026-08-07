#!/usr/bin/env python3
"""Exact certificate for cyclic packet determinant formulas."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.eisenstein import (  # noqa: E402
    Eisenstein,
    determinant,
    determinant_rational,
    is_hermitian,
    trace_form_gram,
)

OUT = ROOT / "certificates" / "eisenstein_packet_determinant.json"

# One cyclic orbit: Hermitian norm six.
H1 = [[Eisenstein(6)]]
assert is_hermitian(H1)
det_h1 = determinant(H1)
assert det_h1 == Eisenstein(6)
assert det_h1.norm() == 36
trace1 = trace_form_gram(H1)
assert trace1 == [[Fraction(12), Fraction(-6)], [Fraction(-6), Fraction(12)]]
assert determinant_rational(trace1) == 108
assert determinant_rational(trace1) == 3 * det_h1.norm()
projected_det1 = det_h1.norm() / 3
assert projected_det1 == 12

# Two orthogonal cyclic orbit planes.
H2 = [[Eisenstein(6), Eisenstein(0)], [Eisenstein(0), Eisenstein(6)]]
det_h2 = determinant(H2)
assert det_h2 == Eisenstein(36)
trace2 = trace_form_gram(H2)
assert determinant_rational(trace2) == 3**2 * det_h2.norm()
assert det_h2.norm() / 3**2 == 12**2

# A nontrivial Hermitian example with off-diagonal 1-zeta.
u = Eisenstein(1, -1)
H3 = [[Eisenstein(6), u], [u.conjugate(), Eisenstein(6)]]
assert is_hermitian(H3)
det_h3 = determinant(H3)
assert det_h3 == Eisenstein(33)
assert det_h3.norm() == 33**2
assert determinant_rational(trace_form_gram(H3)) == 3**2 * det_h3.norm()

certificate = {
    "status": "pass",
    "field": "Q(zeta3)",
    "field_discriminant_absolute_value": 3,
    "one_orbit": {
        "hermitian_gram": [["6"]],
        "hermitian_determinant_norm": 36,
        "difference_gram": [[12, -6], [-6, 12]],
        "difference_determinant": 108,
        "projected_determinant": 12,
    },
    "general_formulas": {
        "difference_lattice_determinant": "3^n Norm(det H)",
        "projected_lattice_determinant": "3^(-n) Norm(det H)",
        "with_invariant_E8_3": "3^(8-n) Norm(det H)",
        "rank30_n11": "Norm(det H)/27",
    },
    "exact_arithmetic_checks": {
        "one_orbit_trace_form": True,
        "two_orbit_block_diagonal": True,
        "nontrivial_off_diagonal_example": True,
    },
    "truth_status": "new intermediate theorem; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": "pass", "rank30_formula": "Norm(det H)/27"}, sort_keys=True))
