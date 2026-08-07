#!/usr/bin/env python3
"""Exact finite-field denominator sieve for low-height twist sections.

The base curve is the certified Kumar--Shioda split-E8 surface

  y^2 = x^3 + t^2 x^2 + a4(t) x + a6(t).

For each quadratic character supported on a pair from
{-B,...,B,infinity}, this script exhausts every natural-degree polynomial
section of

  d(t) y(t)^2 = x(t)^3 + t^2 x(t)^2 + a4(t)x(t)+a6(t)

with deg x <= 2 and deg y <= 2 over F_p.  If no such section exists modulo a
good prime p, then any rational section in this ansatz must have at least one
coefficient with p in its denominator.  Multiple obstruction primes therefore
force a large common denominator.  This is a visibility theorem for the stated
ansatz, not a proof that the full twist Mordell--Weil group is zero.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path

MU = 9699690
RATIONALS = {
    "p2": Fraction(146156773903879871001810589, 2**9 * 3 * MU**2),
    "p1": -Fraction(24909805041567866985469379779685360019313, 2**20 * MU**3),
    "p0": Fraction(14921071761102637668643191215755039801471771138867387, 2**23 * 3 * MU**4),
    "q4": -Fraction(2243374456559366834339, 2**5 * MU**2),
    "q3": Fraction(430800343129403388346226518246078567, 2**11 * MU**3),
    "q2": Fraction(72555101947649011127391733034984158462573146409905769, 2**22 * 3**2 * MU**4),
    "q1": -Fraction(1288109930551729133820743237846836849158406377255698116491924530489, 2**29 * 3 * MU**5),
    "q0": Fraction(8827176793323619929427303381485459401911918837196838709750423283443360357992650203, 2**42 * 3**3 * MU**6),
}


def mod_fraction(value: Fraction, p: int) -> int:
    if value.denominator % p == 0:
        raise ZeroDivisionError(f"prime {p} divides model denominator")
    return (value.numerator % p) * pow(value.denominator % p, -1, p) % p


def model_coefficients(p: int) -> tuple[int, ...]:
    return tuple(mod_fraction(RATIONALS[k], p) for k in
                 ("p0", "p1", "p2", "q0", "q1", "q2", "q3", "q4"))


def rhs_coefficients(A: int, B: int, C: int, p: int, coeffs: tuple[int, ...]) -> tuple[int, ...]:
    """Coefficients in ascending t-order for x^3+t^2*x^2+a4*x+a6."""
    p0, p1, p2, q0, q1, q2, q3, q4 = coeffs
    r = [0] * 7
    r[0] = C**3
    r[1] = 3*C*C*B
    r[2] = 3*C*C*A + 3*C*B*B
    r[3] = B**3 + 6*A*B*C
    r[4] = 3*A*B*B + 3*C*A*A
    r[5] = 3*B*A*A
    r[6] = A**3
    r[2] += C*C
    r[3] += 2*C*B
    r[4] += B*B + 2*C*A
    r[5] += 2*B*A
    r[6] += A*A
    r[0] += p0*C
    r[1] += p0*B + p1*C
    r[2] += p0*A + p1*B + p2*C
    r[3] += p1*A + p2*B
    r[4] += p2*A
    r[0] += q0
    r[1] += q1
    r[2] += q2
    r[3] += q3
    r[4] += q4
    r[5] += 1
    return tuple(x % p for x in r)


def evaluate(poly: tuple[int, ...], x: int, p: int) -> int:
    out = 0
    for c in reversed(poly):
        out = (out*x + c) % p
    return out


def divide_by_linear(poly: tuple[int, ...], root: int, p: int) -> tuple[int, ...]:
    n = len(poly) - 1
    q = [0] * n
    q[-1] = poly[-1] % p
    for k in range(n - 2, -1, -1):
        q[k] = (poly[k + 1] + root*q[k + 1]) % p
    assert (poly[0] + root*q[0]) % p == 0
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return tuple(q)


def square_roots(value: int, p: int) -> list[int]:
    return [x for x in range(p) if x*x % p == value % p]


def polynomial_square_roots(poly: tuple[int, ...], p: int) -> list[tuple[int, ...]]:
    q = list(poly)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    if q == [0]:
        return [(0,)]
    degree = len(q) - 1
    if degree % 2:
        return []
    if degree == 0:
        return [(f,) for f in square_roots(q[0], p)]
    if degree == 2:
        q0, q1, q2 = q
        out = []
        for d in square_roots(q2, p):
            if d == 0:
                continue
            f = q1 * pow(2*d % p, -1, p) % p
            if f*f % p == q0:
                out.append((f, d))
        return out
    if degree == 4:
        q0, q1, q2, q3, q4 = q
        out = []
        for d in square_roots(q4, p):
            if d == 0:
                continue
            e = q3 * pow(2*d % p, -1, p) % p
            f = (q2 - e*e) * pow(2*d % p, -1, p) % p
            if 2*e*f % p == q1 and f*f % p == q0:
                out.append((f, e, d))
        return out
    return []


def pair_label(a: int | None, b: int | None) -> str:
    def one(x): return "infinity" if x is None else str(x)
    return f"{one(a)}|{one(b)}"


def scan_prime(p: int, bound: int) -> dict:
    coeffs = model_coefficients(p)
    rhs_table = [
        ((A, B, C), rhs_coefficients(A, B, C, p, coeffs))
        for A in range(p) for B in range(p) for C in range(p)
    ]
    supports: list[int | None] = list(range(-bound, bound + 1)) + [None]
    if len({x % p for x in supports if x is not None}) != len(supports) - 1:
        raise ValueError(f"support collision modulo {p}")
    pairs = []
    total_sections = 0
    for a, b in combinations(supports, 2):
        aa = None if a is None else a % p
        bb = None if b is None else b % p
        solutions = []
        for abc, rhs in rhs_table:
            if aa is None or bb is None:
                root = bb if aa is None else aa
                if evaluate(rhs, root, p):
                    continue
                quotient = divide_by_linear(rhs, root, p)
            else:
                if evaluate(rhs, aa, p) or evaluate(rhs, bb, p):
                    continue
                quotient = divide_by_linear(divide_by_linear(rhs, aa, p), bb, p)
            roots = polynomial_square_roots(quotient, p)
            if roots:
                solutions.append({"x_coefficients_ABC": list(abc),
                                  "y_coefficients_ascending": [list(r) for r in roots]})
        total_sections += len(solutions)
        pairs.append({"pair": pair_label(a, b),
                      "x_class_count": len(solutions),
                      "samples": solutions[:4]})
    nonempty = [row for row in pairs if row["x_class_count"]]
    return {
        "prime": p,
        "bound": bound,
        "pair_count": len(pairs),
        "nonempty_pair_count": len(nonempty),
        "total_x_classes": total_sections,
        "nonempty_pairs": nonempty,
        "all_pairs_empty": not nonempty,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="23,29,31,37,41,47,53,61,71,73,79")
    parser.add_argument("--bound", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primes = [int(x) for x in args.primes.split(",") if x.strip()]
    runs = [scan_prime(p, args.bound) for p in primes]
    universal_obstructions = [row["prime"] for row in runs if row["all_pairs_empty"]]
    out = {
        "status": "pass",
        "ansatz": {
            "twist_equation": "d(t)y(t)^2=x(t)^3+t^2x(t)^2+a4(t)x(t)+a6(t)",
            "x_degree_bound": 2,
            "y_degree_bound": 2,
            "branch_support": [str(x) for x in range(-args.bound, args.bound + 1)] + ["infinity"],
        },
        "runs": runs,
        "universal_obstruction_primes": universal_obstructions,
        "forced_common_denominator_divisor": prod(universal_obstructions),
        "theorem": (
            "For every branch pair in the stated support, any Q-rational section in the stated "
            "polynomial ansatz must have common coefficient denominator divisible by every "
            "universal obstruction prime."
        ),
        "truth_note": (
            "This is an exact denominator/visibility obstruction for one low-height ansatz. "
            "It is not a Mordell-Weil rank upper bound and does not exclude rational-function sections."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": out["status"],
        "universal_obstruction_primes": universal_obstructions,
        "forced_common_denominator_divisor": out["forced_common_denominator_divisor"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
