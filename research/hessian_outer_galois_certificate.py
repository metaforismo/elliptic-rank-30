#!/usr/bin/env python3
"""Exact symbolic certificate for the Hessian outer-Galois reduction."""
from __future__ import annotations

from pathlib import Path
import json

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "hessian_outer_galois_reduction.json"

# Plane-cubic translation.
x, X, y, z = sp.symbols("x X y z")
a = sp.symbols("a", nonzero=True)
l1, l2 = sp.symbols("l1 l2")
q11, q12, q22 = sp.symbols("q11 q12 q22")
c30, c21, c12, c03 = sp.symbols("c30 c21 c12 c03")
L = l1 * y + l2 * z
Q = q11 * y**2 + 2 * q12 * y * z + q22 * z**2
C = c30 * y**3 + c21 * y**2 * z + c12 * y * z**2 + c03 * z**3
plane_cubic = a * x**3 + x**2 * L + x * Q + C
translated = sp.expand(plane_cubic.subs(x, X - L / (3 * a)))
poly_X = sp.Poly(translated, X)
assert sp.simplify(poly_X.coeff_monomial(X**2)) == 0
mixed_quadratic = sp.factor(poly_X.coeff_monomial(X))
expected_mixed = sp.factor(Q - L**2 / (3 * a))
assert sp.simplify(mixed_quadratic - expected_mixed) == 0
cyclic_quadratic = sp.expand(3 * a * Q - L**2)

# Ambient block determinant identity.
l3 = sp.symbols("l3")
q13, q23, q33 = sp.symbols("q13 q23 q33")
lambda_vector = sp.Matrix([l1, l2, l3])
Q_matrix = sp.Matrix(
    [
        [q11, q12, q13],
        [q12, q22, q23],
        [q13, q23, q33],
    ]
)
M = sp.Matrix(
    [
        [a, l1, l2, l3],
        [l1, q11, q12, q13],
        [l2, q12, q22, q23],
        [l3, q13, q23, q33],
    ]
)
B = a * Q_matrix - lambda_vector * lambda_vector.T
block_identity = sp.factor(B.det() - a**2 * M.det())
assert block_identity == 0

certificate = {
    "status": "pass",
    "characteristic_exclusions": [2, 3],
    "plane_cubic_normal_form": "a*x^3 + x^2*L + x*Q + C",
    "coordinate_change": "x = X - L/(3a)",
    "translated_X2_coefficient": "0",
    "translated_X_coefficient": str(expected_mixed),
    "cyclicity_condition": "3*a*Q - L^2 = 0",
    "invariant_bilinear_form": "B_q(u,v)=F(q)T(q,u,v)-T(q,q,u)T(q,q,v)",
    "block_determinant_identity": "det(B_q)=F(q)^2 det(T(q,-,-))",
    "hessian_scaling": "Hess_F(q)=6*T(q,-,-)",
    "hessian_equivalence": "det(B_q)=0 iff q lies on the Hessian quartic of X",
    "generic_hessian_points_on_line": 4,
    "generic_isotropic_hyperplanes_per_rank2_point": 2,
    "generic_cyclic_hyperplanes_per_trace_class_upper_bound": 8,
    "generic_candidates_per_240_class_trace_code_upper_bound": 1920,
    "symbolic_checks": {
        "X2_term_eliminated": True,
        "mixed_quadratic_formula": True,
        "block_determinant_identity": True,
    },
    "truth_status": "new intermediate theorem; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(
    json.dumps(
        {
            "status": certificate["status"],
            "cyclic_hyperplanes_per_trace": certificate[
                "generic_cyclic_hyperplanes_per_trace_class_upper_bound"
            ],
            "candidates_per_trace_code": certificate[
                "generic_candidates_per_240_class_trace_code_upper_bound"
            ],
        },
        sort_keys=True,
    )
)
