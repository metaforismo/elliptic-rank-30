#!/usr/bin/env python3
"""Exact certificate for a genus-zero three-channel quadratic rank packet.

The family is

    E_t: y^2 = x^3 - (t^2 + 1)/(t + 1) * x + t.

Let d1=t and d2=2t/(t+1). The twists d1*y^2=f_t(x),
d2*y^2=f_t(x), and d1*d2*y^2=f_t(x) contain points with
x-coordinates 0, 1, and -1, respectively. Their squareclasses have total
geometric branch support {0,-1,infinity}; the biquadratic base therefore has
genus zero. The three points occupy the three distinct non-trivial Galois
characters.

This module verifies every polynomial identity with exact integer arithmetic,
checks the branch-code calculation, verifies a rational parametrisation of the
base, and proves that the three generic twist sections are non-torsion by a
Lutz--Nagell specialization at t=-3.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

Poly = Tuple[int, ...]
Point = Optional[Tuple[Fraction, Fraction]]


def trim(a: Iterable[int]) -> Poly:
    values = list(a)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [0])


def padd(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return trim((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n))


def pneg(a: Poly) -> Poly:
    return tuple(-x for x in a)


def psub(a: Poly, b: Poly) -> Poly:
    return padd(a, pneg(b))


def pmul(a: Poly, b: Poly) -> Poly:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def pscale(a: Poly, c: int) -> Poly:
    return trim(c * x for x in a)


def ppow(a: Poly, n: int) -> Poly:
    if n < 0:
        raise ValueError("negative polynomial exponent")
    out: Poly = (1,)
    base = a
    k = n
    while k:
        if k & 1:
            out = pmul(out, base)
        base = pmul(base, base)
        k >>= 1
    return out


def rank_f2(rows: Sequence[Sequence[int]]) -> int:
    if not rows:
        return 0
    matrix = [[entry & 1 for entry in row] for row in rows]
    nrows, ncols = len(matrix), len(matrix[0])
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, nrows) if matrix[r][col]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for r in range(nrows):
            if r != rank and matrix[r][col]:
                matrix[r] = [x ^ y for x, y in zip(matrix[r], matrix[rank])]
        rank += 1
    return rank


def discriminant_short(a: int, b: int) -> int:
    return -16 * (4 * a**3 + 27 * b**2)


def ec_add(a: int, p: Point, q: Point) -> Point:
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2 and y1 == -y2:
        return None
    if p == q:
        if y1 == 0:
            return None
        slope = (3 * x1 * x1 + a) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return (x3, y3)


def ec_mul(a: int, n: int, p: Point) -> Point:
    if n < 0:
        if p is None:
            return None
        return ec_mul(a, -n, (p[0], -p[1]))
    out: Point = None
    base = p
    k = n
    while k:
        if k & 1:
            out = ec_add(a, out, base)
        base = ec_add(a, base, base)
        k >>= 1
    return out


@dataclass(frozen=True)
class NontorsionWitness:
    twist: int
    curve_a: int
    curve_b: int
    point: Tuple[int, int]
    multiple: int
    multiple_x: Fraction


def verify_packet() -> dict:
    t: Poly = (0, 1)
    t_plus_one: Poly = (1, 1)
    t2_plus_one: Poly = (1, 0, 1)

    def rhs_numerator(c: int) -> Poly:
        return psub(pmul(padd((c**3,), t), t_plus_one), pscale(t2_plus_one, c))

    rhs0 = rhs_numerator(0)
    rhs1 = rhs_numerator(1)
    rhsm1 = rhs_numerator(-1)
    assert rhs0 == pmul(t, t_plus_one)
    assert rhs1 == pscale(t, 2)
    assert rhsm1 == pscale(ppow(t, 2), 2)

    d1_vector = (1, 0, 1)
    d2_vector = (1, 1, 0)
    d3_vector = tuple(x ^ y for x, y in zip(d1_vector, d2_vector))
    assert d3_vector == (0, 1, 1)
    character_rank = rank_f2([d1_vector, d2_vector])
    assert character_rank == 2
    branch_support_size = 3
    genus = 1 + 2 ** (character_rank - 2) * (branch_support_size - 4)
    assert genus == 0

    den: Poly = (1, 0, 1)
    v_num: Poly = (-1, -2, 1)
    w_num: Poly = (1, -2, -1)
    assert padd(ppow(v_num, 2), ppow(w_num, 2)) == pscale(ppow(den, 2), 2)
    lhs = pmul(ppow(v_num, 2), padd(ppow(w_num, 2), ppow(v_num, 2)))
    rhs = pscale(pmul(ppow(den, 2), ppow(v_num, 2)), 2)
    assert lhs == rhs

    witnesses = []
    for d, x_on_d_model, multiple in [(-3, 0, 2), (3, 1, 3), (-9, -1, 3)]:
        a = d * d * 5
        b = d * d * d * (-3)
        assert discriminant_short(a, b) != 0
        point = (Fraction(d * x_on_d_model), Fraction(d * d))
        assert point[1] * point[1] == point[0] ** 3 + a * point[0] + b
        multiple_point = ec_mul(a, multiple, point)
        assert multiple_point is not None
        assert multiple_point[0].denominator != 1
        witnesses.append(
            NontorsionWitness(
                twist=d,
                curve_a=a,
                curve_b=b,
                point=(int(point[0]), int(point[1])),
                multiple=multiple,
                multiple_x=multiple_point[0],
            )
        )

    assert witnesses[0].multiple_x == Fraction(25, 4)
    assert witnesses[1].multiple_x == Fraction(1479, 49)
    assert witnesses[2].multiple_x == Fraction(13077, 121)

    return {
        "status": "pass",
        "theorem": "genus-zero three-channel quadratic rank packet",
        "family": {
            "base_field": "Q(t)",
            "equation": "y^2 = x^3 - (t^2+1)/(t+1) * x + t",
            "excluded_parameters": ["t=-1", "zeros of the discriminant"],
        },
        "characters": {
            "d1": "t",
            "d2": "2*t/(t+1)",
            "d1d2": "2*t^2/(t+1)",
            "geometric_dimension": character_rank,
            "branch_support": ["0", "-1", "infinity"],
            "branch_support_size": branch_support_size,
            "cover_genus": genus,
        },
        "twist_sections": [
            {"character": "d1", "x": "0", "twist_y": "1"},
            {"character": "d2", "x": "1", "twist_y": "1"},
            {"character": "d1d2", "x": "-1", "twist_y": "1"},
        ],
        "base_change_parametrisation": {
            "v": "(r^2-2*r-1)/(r^2+1)",
            "w": "(1-2*r-r^2)/(r^2+1)",
            "u": "v/w",
            "t": "u^2",
            "sqrt_d1": "u",
            "sqrt_d2": "v",
            "identity": "v^2+w^2=2",
        },
        "nontorsion_specialization": {
            "parameter": "t=-3",
            "base_curve": "y^2=x^3+5*x-3",
            "witnesses": [
                {
                    "twist": w.twist,
                    "curve": f"y^2=x^3+({w.curve_a})*x+({w.curve_b})",
                    "point": list(w.point),
                    "multiple": w.multiple,
                    "multiple_x": f"{w.multiple_x.numerator}/{w.multiple_x.denominator}",
                    "lutz_nagell_conclusion": "non-torsion",
                }
                for w in witnesses
            ],
        },
        "rank_conclusion": "rank E(L) >= rank E(Q(t)) + 3",
        "independence_reason": "the three non-torsion sections lie in distinct non-trivial characters of Gal(L/Q(t))",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_packet()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
