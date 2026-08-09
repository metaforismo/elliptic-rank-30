#!/usr/bin/env python3
"""Certify trisection-packet rank directly on specialized plane cubics.

This is independent of any global cubic-to-Weierstrass transformation.  It
reduces the plane cubic and all displayed rational points modulo good primes,
forms the packet differences with the chord-tangent group law, and accumulates
exact ranks in E(F_p)/ell E(F_p).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

from plane_cubic_finite_reduction import (
    PlaneCubicGroup,
    Point,
    rational_projective_mod,
    row_rank,
)

BASE_POINTS: list[list[Fraction]] = [
    [Fraction(-2), Fraction(-1), Fraction(1)],
    [Fraction(-4), Fraction(2), Fraction(1)],
    [Fraction(1), Fraction(-3), Fraction(1)],
    [Fraction(4), Fraction(-3), Fraction(1)],
    [Fraction(2), Fraction(-2), Fraction(1)],
    [Fraction(0), Fraction(3), Fraction(1)],
    [Fraction(3), Fraction(4), Fraction(1)],
    [Fraction(-3), Fraction(0), Fraction(1)],
    [Fraction(5264, 13547), Fraction(-35295, 13547), Fraction(1)],
]
ORIGIN = BASE_POINTS[8]

KNOWN_T = Fraction(1427161, 111293)
KNOWN_ROOTS = [
    Fraction(-16253, 2248),
    Fraction(-7838, 67937),
    Fraction(34, 19),
]

POINT_FORMULAS = {
    "C1": {
        "den": [-118794, 154067, 179873, 437604],
        "x": [72912, 289296, -2367708],
        "y": [-356382, 535113, 828915, -1054896],
    },
    "C2": {
        "den": [-17044, 6464, -23081, 13816],
        "x": [46142, 104745, -91352],
        "y": [-51132, 65534, 35502, -49904],
    },
    "Q": {
        "den": [5590868988, -13577828351, 8426176988],
        "x": [-12673653216, 25187117592, -7840260876],
        "y": [-4236741228, 15734911491, -15419156388],
    },
}


def qstr(value: Fraction | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def as_fraction(value) -> Fraction:
    """Coerce Python or Sage rational values to ``fractions.Fraction``."""
    if isinstance(value, Fraction):
        return value
    numerator = value.numerator() if callable(getattr(value, "numerator", None)) else value.numerator
    denominator = value.denominator() if callable(getattr(value, "denominator", None)) else value.denominator
    return Fraction(int(numerator), int(denominator))


def evaluate(coefficients: Sequence[int], s: Fraction) -> Fraction:
    s = as_fraction(s)
    return sum((Fraction(c) * s**i for i, c in enumerate(coefficients)), Fraction(0))


def plane_point(s: Fraction, kind: str) -> list[Fraction]:
    s = as_fraction(s)
    formula = POINT_FORMULAS[kind]
    den = evaluate(formula["den"], s)
    if den == 0:
        raise ZeroDivisionError(f"{kind} parameterization pole at {s}")
    return [evaluate(formula["x"], s) / den, evaluate(formula["y"], s) / den, Fraction(1)]


def packet_points_mod(group: PlaneCubicGroup, roots: Sequence[Fraction]) -> tuple[list[Point], list[str], list[dict]]:
    points = [rational_projective_mod(P, group.p) for P in BASE_POINTS[:8]]
    labels = [f"B{i}" for i in range(8)]
    packet_records = []
    for j, s in enumerate(roots):
        raw = {kind: plane_point(s, kind) for kind in ("C1", "C2", "Q")}
        reduced = {kind: rational_projective_mod(P, group.p) for kind, P in raw.items()}
        for kind, P in reduced.items():
            if P not in group.point_set:
                raise ValueError(("packet point off reduced cubic", j, kind, group.p, P))
        D1 = group.add(reduced["C1"], group.neg(reduced["Q"]))
        D2 = group.add(reduced["C2"], group.neg(reduced["Q"]))
        points.extend([D1, D2])
        labels.extend([f"D1_{j}", f"D2_{j}"])
        packet_records.append({
            "root": qstr(s),
            "plane_points": {kind: [qstr(x) for x in P] for kind, P in raw.items()},
        })
    return points, labels, packet_records


def certify_packet(
    *,
    t: Fraction,
    roots: Sequence[Fraction],
    ell: int,
    primes: Iterable[int],
    target_rank: int | None = None,
    expected_orders: dict[int, int] | None = None,
) -> dict:
    target_rank = len(BASE_POINTS[:8]) + 2 * len(roots) if target_rank is None else target_rank
    rows = [[] for _ in range(len(BASE_POINTS[:8]) + 2 * len(roots))]
    records = []
    torsion_witness = None
    labels = None
    packet_records = None
    failures = []

    for p in primes:
        try:
            group = PlaneCubicGroup.create(p, t.numerator, t.denominator, ORIGIN)
            if expected_orders and p in expected_orders and len(group.points) != expected_orders[p]:
                raise AssertionError(("group order mismatch", p, len(group.points), expected_orders[p]))
            reduced, labels, packet_records = packet_points_mod(group, roots)
            vectors, dimension = group.quotient_vectors_points(reduced, ell)
        except (ValueError, ZeroDivisionError, AssertionError) as exc:
            failures.append({"prime": p, "error": repr(exc)})
            continue

        order = len(group.points)
        if torsion_witness is None and order % ell:
            torsion_witness = {"prime": p, "group_order": order}
        if dimension:
            records.append({"prime": p, "group_order": order, "quotient_dimension": dimension})
            for i, vector in enumerate(vectors):
                rows[i].extend(vector)
        rank = row_rank(rows, ell)
        if rank >= target_rank and torsion_witness is not None:
            break

    rank = row_rank(rows, ell)
    return {
        "ell": ell,
        "rank": rank,
        "target_rank": target_rank,
        "rows": rows,
        "primes": records,
        "torsion_witness": torsion_witness,
        "labels": labels,
        "packet_records": packet_records,
        "failures": failures,
    }


def canonical_hash(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def load_tangent_certificate(path: Path) -> tuple[Fraction, list[Fraction], str]:
    data = json.loads(path.read_text())
    t = Fraction(data["common_fibre_parameter"])
    roots = [Fraction(x) for x in data["three_source_roots"]]
    return t, roots, data["record_sha256"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("known", "tangent"), required=True)
    parser.add_argument("--ell", type=int, default=2)
    parser.add_argument("--primes", default="")
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.mode == "known":
        t, roots = KNOWN_T, KNOWN_ROOTS
        source = {"kind": "known split fibre", "t": qstr(t)}
        if args.primes:
            primes = [int(x) for x in args.primes.split(",") if x]
        else:
            primes = [43, 53, 61, 67, 101, 191]
        expected_orders = {43: 52, 53: 60, 61: 76, 67: 73, 101: 104, 191: 192}
    else:
        if args.certificate is None:
            raise SystemExit("--certificate is required in tangent mode")
        t, roots, parent_hash = load_tangent_certificate(args.certificate)
        source = {"kind": "tangent-derived split fibre", "parent_certificate_sha256": parent_hash, "t": qstr(t)}
        if args.primes:
            primes = [int(x) for x in args.primes.split(",") if x]
        else:
            primes = [43, 53, 61, 67, 101, 191]
        expected_orders = None

    result = certify_packet(
        t=t,
        roots=roots,
        ell=args.ell,
        primes=primes,
        target_rank=14,
        expected_orders=expected_orders,
    )
    payload = {
        "schema_version": 1,
        "certificate_id": f"small-split-e8-three-trisection-plane-cubic-{args.mode}-ell{args.ell}-v1",
        "claim_status": "exact finite-reduction diagnostic",
        "source": source,
        "t": qstr(t),
        "roots": [qstr(s) for s in roots],
        "point_count": 14,
        "result": result,
        "truth_note": (
            "Full row rank together with the stored torsion witness proves independence modulo torsion. "
            "A smaller rank is only a lower bound/diagnostic and does not prove an upper bound."
        ),
    }
    payload["record_sha256"] = canonical_hash(payload)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
