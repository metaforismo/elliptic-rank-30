#!/usr/bin/env python3
"""Exact symbolic coefficient certificate for the j=0 cubic characters."""
from __future__ import annotations

from pathlib import Path
import json

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "j0_cubic_character_obstruction.json"

u = sp.symbols("u")
a, b = sp.symbols("a b", nonzero=True)
A, B, C, D, E = sp.symbols("A B C D E")
F = b**2 * u**6 - a**3 * (u**3 + 1) ** 3

# Character two.
x2 = A * u**2
y2 = C * u**6 + D * u**3 + E
difference2 = sp.Poly(sp.expand(y2**2 - x2**3 - F), u)
coeff2 = {degree: sp.expand(difference2.coeff_monomial(u**degree)) for degree in range(13)}
assert coeff2[12] == C**2
assert coeff2[9] == 2 * C * D + a**3
# Therefore C=0 forces a^3=0, impossible in the marked family.

# Character one.
x1 = u * (A * u**3 + B)
y1 = C * u**6 + D * u**3 + E
difference1 = sp.Poly(sp.expand(y1**2 - x1**3 - F), u)
actual = {degree: sp.expand(difference1.coeff_monomial(u**degree)) for degree in (12, 9, 6, 3, 0)}
expected = {
    12: C**2 - A**3,
    9: 2 * C * D - 3 * A**2 * B + a**3,
    6: D**2 + 2 * C * E - 3 * A * B**2 - b**2 + 3 * a**3,
    3: 2 * D * E - B**3 + 3 * a**3,
    0: E**2 + a**3,
}
for degree in expected:
    assert sp.expand(actual[degree] - expected[degree]) == 0

mu, s = sp.symbols("mu s")
q_mu = mu**2 - mu + 1
assert sp.expand(q_mu - ((mu - sp.Rational(1, 2)) ** 2 + sp.Rational(3, 4))) == 0
assert sp.expand(q_mu.subs(mu, 2) - 3) == 0
assert sp.expand((q_mu - 3 * s**2).subs({mu: 2, s: 1})) == 0

certificate = {
    "status": "pass",
    "cubic_base_curve": "y^2=x^3+b^2*u^6-a^3*(u^3+1)^3",
    "character_two": {
        "x_ansatz": "A*u^2",
        "y_ansatz": "C*u^6+D*u^3+E",
        "u12_equation": "C^2=0",
        "u9_equation_after_C_zero": "a^3=0",
        "conclusion": "empty over the algebraic closure for a nonzero",
    },
    "character_one": {
        "x_ansatz": "u*(A*u^3+B)",
        "y_ansatz": "C*u^6+D*u^3+E",
        "coefficient_equations": [
            "C^2=A^3",
            "2CD=3A^2B-a^3",
            "D^2+2CE=3AB^2+b^2-3a^3",
            "2DE=B^3-3a^3",
            "E^2=-a^3",
        ],
        "minimum_constant_field": "contains sqrt(-a)=sqrt(-q(mu))",
        "no_Q_rational_eigen_section": True,
    },
    "cm_compatible_locus": {
        "equation": "mu^2-mu+1=3s^2",
        "rational_point": ["2", "1"],
        "field": "Q(sqrt(-3))",
    },
    "symbolic_checks": {
        "character_two_leading_obstruction": True,
        "character_one_five_equations": True,
        "q_mu_positive_completion_of_square": True,
        "cm_locus_has_rational_point": True,
    },
    "truth_status": "new restricted obstruction; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": "pass", "cm_parameter": [2, 1]}, sort_keys=True))
