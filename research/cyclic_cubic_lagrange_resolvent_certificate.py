#!/usr/bin/env python3
"""Exact symbolic certificate for the cyclic cubic Lagrange resolvent."""
from __future__ import annotations

from pathlib import Path
import json

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "cyclic_cubic_lagrange_resolvent.json"

A, B, C, X = sp.symbols("A B C X")
p = B - A**2 / 3
q = 2 * A**3 / 27 - A * B / 3 + C
original = (X - A / 3) ** 3 + A * (X - A / 3) ** 2 + B * (X - A / 3) + C
assert sp.expand(original - (X**3 + p * X + q)) == 0

P, Q, S = sp.symbols("p q s")
delta = -4 * P**3 - 27 * Q**2
sqrt_minus_three = sp.symbols("sqrt_minus_three")
d_plus = (-27 * Q + 3 * S * sqrt_minus_three) / 2
d_minus = (-27 * Q - 3 * S * sqrt_minus_three) / 2

# Work modulo sqrt(-3)^2=-3 and s^2=Delta.
product = sp.expand(d_plus * d_minus).subs(sqrt_minus_three**2, -3)
product = sp.expand(product).subs(S**2, delta)
assert sp.factor(product + 27 * P**3) == 0
assert sp.factor(d_plus + d_minus + 27 * Q) == 0

resolvent_discriminant = sp.expand((27 * Q) ** 2 + 108 * P**3)
assert sp.factor(resolvent_discriminant + 27 * delta) == 0

certificate = {
    "status": "pass",
    "characteristic_exclusions": [2, 3],
    "depressed_cubic": {
        "p": "B-A^2/3",
        "q": "2A^3/27-AB/3+C",
        "discriminant": "-4p^3-27q^2",
    },
    "lagrange_cubes": {
        "R_cubed": "(-27q+3s*sqrt(-3))/2",
        "S_cubed": "(-27q-3s*sqrt(-3))/2",
        "sum": "-27q",
        "product": "-27p^3",
    },
    "kummer_hash": "unordered pair {[d],[d]^-1} in K(zeta3)^*/K(zeta3)^{*3}",
    "symbolic_checks": {
        "depression_formula": True,
        "resolvent_sum": True,
        "resolvent_product_mod_relations": True,
        "resolvent_discriminant_equals_minus_27_Delta": True,
    },
    "truth_status": "exact algebraic reduction; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": "pass", "checks": certificate["symbolic_checks"]}, sort_keys=True))
