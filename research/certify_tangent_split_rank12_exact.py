#!/usr/bin/env python3
"""Build an exact rank-12 packet certificate on the second split fibre.

The proof has two parts:
  * two exact trace relations give subgroup rank at most 12;
  * exact finite reductions modulo ell=5 give rank at least 12.

All rational packet differences are computed directly on the plane cubic with
its specified rational origin, avoiding any giant global Weierstrass map.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from certify_three_trisection_plane_cubic import (
    BASE_POINTS,
    ORIGIN,
    plane_point,
    load_tangent_certificate,
    certify_packet,
)
from plane_cubic_finite_reduction import F_COEFFS, G_COEFFS, row_rank

QPoint = tuple[Fraction, Fraction, Fraction]


def qstr(q0: Fraction | int) -> str:
    q = Fraction(q0)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def qnormalize(values: Sequence[Fraction | int]) -> QPoint:
    v = [Fraction(x) for x in values]
    pivot = next((v[i] for i in range(len(v)-1, -1, -1) if v[i]), None)
    if pivot is None:
        raise ValueError("zero projective vector")
    return tuple(x / pivot for x in v)  # type: ignore[return-value]


def combined_coefficients(t: Fraction) -> dict[tuple[int, int, int], int]:
    keys = set(F_COEFFS) | set(G_COEFFS)
    return {
        m: t.denominator * F_COEFFS.get(m, 0) + t.numerator * G_COEFFS.get(m, 0)
        for m in keys
    }


def qpoly(coeffs: dict[tuple[int, int, int], int], P: QPoint) -> Fraction:
    x, y, z = P
    return sum(Fraction(c) * x**i * y**j * z**k for (i, j, k), c in coeffs.items())


def qgradient(coeffs: dict[tuple[int, int, int], int], P: QPoint) -> tuple[Fraction, Fraction, Fraction]:
    values = [P[0], P[1], P[2]]
    out = []
    for axis in range(3):
        total = Fraction(0)
        for exps, c in coeffs.items():
            e = exps[axis]
            if e == 0:
                continue
            term = Fraction(c * e)
            for q, exponent in enumerate(exps):
                term *= values[q] ** (exponent - (1 if q == axis else 0))
            total += term
        out.append(total)
    return tuple(out)  # type: ignore[return-value]


class RationalPlaneCubicGroup:
    def __init__(self, coeffs: dict[tuple[int, int, int], int], origin: Sequence[Fraction | int]):
        self.coeffs = coeffs
        self.origin = qnormalize(origin)
        if qpoly(coeffs, self.origin) != 0:
            raise ValueError("origin is off the cubic")

    @staticmethod
    def _add_vectors(P: QPoint, Q: QPoint, sign: int = 1) -> QPoint:
        return tuple(a + sign*b for a, b in zip(P, Q))  # type: ignore[return-value]

    def _tangent_companion(self, P: QPoint) -> QPoint:
        a, b, c = qgradient(self.coeffs, P)
        candidates = []
        if a or b:
            candidates.append((b, -a, Fraction(0)))
        if a or c:
            candidates.append((c, Fraction(0), -a))
        if b or c:
            candidates.append((Fraction(0), c, -b))
        for V in candidates:
            if V != (0, 0, 0) and qnormalize(V) != P:
                return V
        raise AssertionError(("no tangent companion", P, (a, b, c)))

    def third(self, P: QPoint, Q: QPoint) -> QPoint:
        if P != Q:
            plus = qpoly(self.coeffs, self._add_vectors(P, Q, 1))
            minus = qpoly(self.coeffs, self._add_vectors(P, Q, -1))
            A = (plus - minus) / 2
            B = (plus + minus) / 2
            if A == 0 and B == 0:
                raise AssertionError("line is a component")
            R = qnormalize(tuple(B*x - A*y for x, y in zip(P, Q)))
        else:
            V = self._tangent_companion(P)
            plus = qpoly(self.coeffs, self._add_vectors(P, V, 1))
            minus = qpoly(self.coeffs, self._add_vectors(P, V, -1))
            A = (plus + minus) / 2
            B = (plus - minus) / 2
            if B != qpoly(self.coeffs, V):
                raise AssertionError("tangent binary-cubic mismatch")
            R = P if A == 0 else qnormalize(tuple(-B*x + A*y for x, y in zip(P, V)))
        if qpoly(self.coeffs, R) != 0:
            raise AssertionError("third point is off the cubic")
        return R

    def add(self, P: QPoint, Q: QPoint) -> QPoint:
        return self.third(self.third(P, Q), self.origin)

    def neg(self, P: QPoint) -> QPoint:
        return self.third(self.third(self.origin, self.origin), P)

    def mul(self, n: int, P: QPoint) -> QPoint:
        if n < 0:
            return self.mul(-n, self.neg(P))
        R = self.origin
        Q = P
        while n:
            if n & 1:
                R = self.add(R, Q)
            n //= 2
            if n:
                Q = self.add(Q, Q)
        return R

    def linear_combination(self, coefficients: Sequence[int], points: Sequence[QPoint]) -> QPoint:
        if len(coefficients) != len(points):
            raise ValueError("coefficient length mismatch")
        R = self.origin
        for n, P in zip(coefficients, points):
            if n:
                R = self.add(R, self.mul(n, P))
        return R


def independent_indices(rows: Sequence[Sequence[int]], ell: int) -> list[int]:
    chosen: list[int] = []
    current_rank = 0
    for i in range(len(rows)):
        new_rank = row_rank([rows[j] for j in chosen + [i]], ell)
        if new_rank > current_rank:
            chosen.append(i)
            current_rank = new_rank
    return chosen


def canonical_hash(payload: dict) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split-certificate", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--check", type=Path)
    args = ap.parse_args()

    t, roots, split_hash = load_tangent_certificate(args.split_certificate)
    reduction = certify_packet(
        t=t,
        roots=roots,
        ell=5,
        primes=[29, 59, 61, 73, 89, 131, 139, 179, 181, 251, 281, 307,
                311, 313, 317, 331, 337, 349],
        target_rank=12,
    )
    if reduction["ell"] != 5 or reduction["rank"] != 12:
        raise AssertionError("failed to recompute an ell=5 rank-12 reduction certificate")
    if reduction["torsion_witness"] is None:
        raise AssertionError("missing torsion witness")

    coeffs = combined_coefficients(t)
    group = RationalPlaneCubicGroup(coeffs, ORIGIN)
    base = [qnormalize(P) for P in BASE_POINTS[:8]]
    all_points: list[QPoint] = list(base)
    labels = [f"B{i}" for i in range(8)]
    raw_packet = []
    for j, s in enumerate(roots):
        C1 = qnormalize(plane_point(s, "C1"))
        C2 = qnormalize(plane_point(s, "C2"))
        Q = qnormalize(plane_point(s, "Q"))
        for P in (C1, C2, Q):
            if qpoly(coeffs, P) != 0:
                raise AssertionError(("packet point off cubic", j, P))
        D1 = group.add(C1, group.neg(Q))
        D2 = group.add(C2, group.neg(Q))
        all_points.extend([D1, D2])
        labels.extend([f"D1_{j}", f"D2_{j}"])
        raw_packet.append({
            "source_root": qstr(s),
            "C1": [qstr(x) for x in C1],
            "C2": [qstr(x) for x in C2],
            "Q": [qstr(x) for x in Q],
            "D1": [qstr(x) for x in D1],
            "D2": [qstr(x) for x in D2],
        })

    if len(all_points) != 14 or reduction["labels"] != labels:
        raise AssertionError("point-label mismatch")
    if any(qpoly(coeffs, P) != 0 for P in all_points):
        raise AssertionError("displayed point off cubic")

    relation_1 = [-1, -1, -1, -1, -4, 5, 2, -1, 3, 0, 3, 0, 3, 0]
    relation_2 = [-4, 2, 2, -4, -1, 5, -1, 2, 0, 3, 0, 3, 0, 3]
    relation_checks = []
    for vector in (relation_1, relation_2):
        value = group.linear_combination(vector, all_points)
        relation_checks.append(value == group.origin)
    if relation_checks != [True, True]:
        raise AssertionError(("trace relation failure", relation_checks))

    selected = independent_indices(reduction["rows"], 5)
    if len(selected) != 12:
        raise AssertionError(("independent subset size", len(selected), selected))

    payload = {
        "schema_version": 1,
        "certificate_id": "small-split-e8-tangent-fibre-exact-rank12-packet-v1",
        "claim_status": "certified exact subgroup rank 12",
        "curve": {
            "model": "smooth plane cubic with specified rational origin",
            "equation": "sum c_ijk X^i Y^j Z^k = 0",
            "integral_coefficients": {
                f"{i},{j},{k}": str(c) for (i, j, k), c in sorted(coeffs.items())
            },
            "fibre_parameter": qstr(t),
            "origin": [qstr(x) for x in group.origin],
            "nonsingularity_witness": {
                "good_reduction_prime": reduction["torsion_witness"]["prime"],
                "reduced_group_order": reduction["torsion_witness"]["group_order"],
                "reason": "the direct reduction verifier checked that the projective cubic is smooth at every F_p-point",
            },
        },
        "parent_split_certificate_sha256": split_hash,
        "labels": labels,
        "points": {label: [qstr(x) for x in P] for label, P in zip(labels, all_points)},
        "raw_packet_points": raw_packet,
        "exact_point_on_curve_checks": {label: True for label in labels},
        "relations": [
            {
                "name": "trace_relation_C1_minus_Q",
                "coefficients_in_label_order": relation_1,
                "exact_group_law_check": True,
            },
            {
                "name": "trace_relation_C2_minus_Q",
                "coefficients_in_label_order": relation_2,
                "exact_group_law_check": True,
            },
        ],
        "rank_upper_bound": {
            "value": 12,
            "reason": "two Q-linearly independent exact integral relations among 14 points",
        },
        "rank_lower_bound": {
            "value": 12,
            "ell": 5,
            "finite_reduction_rank": reduction["rank"],
            "torsion_witness": reduction["torsion_witness"],
            "used_primes": reduction["primes"],
            "rows": reduction["rows"],
            "reduction_failures": reduction["failures"],
            "selected_independent_indices": selected,
            "selected_independent_labels": [labels[i] for i in selected],
            "reason": (
                "the displayed classes have rank 12 in a product of E(F_p)/5E(F_p); "
                "a good prime with group order prime to 5 excludes rational 5-torsion, "
                "so infinite descent proves Z-independence modulo torsion"
            ),
        },
        "displayed_subgroup_rank": 12,
        "truth_note": (
            "This is an unconditional rank-12 certificate for the displayed subgroup on one explicit fibre. "
            "It is not a rank-30 certificate and does not determine the full Mordell-Weil rank."
        ),
    }
    payload["record_sha256"] = canonical_hash(payload)

    if args.check:
        previous = json.loads(args.check.read_text())
        if previous != payload:
            raise SystemExit("certificate mismatch")
        print("VERIFIED exact rank-12 packet sha256=" + payload["record_sha256"])
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "displayed_subgroup_rank": payload["displayed_subgroup_rank"],
        "selected_independent_labels": payload["rank_lower_bound"]["selected_independent_labels"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
