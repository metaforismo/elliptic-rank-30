#!/usr/bin/env sage-python
"""Independent SageMath replay of the tangent-derived rank-12 packet.

This verifier is deliberately separate from the direct plane-cubic group-law
certificate. Sage constructs a Weierstrass model from the specialized plane
cubic, maps the displayed points, verifies the two trace relations exactly,
and proves that twelve selected points have independent classes modulo 5
using reductions at automatically selected good primes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (  # type: ignore
    EllipticCurve,
    EllipticCurve_from_cubic,
    GF,
    PolynomialRing,
    QQ,
    matrix,
    prime_range,
    version,
)

from certify_three_trisection_plane_cubic import BASE_POINTS, ORIGIN, plane_point


def canonical_hash(payload: dict) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def coords(point) -> list[str]:
    return ["0"] if point.is_zero() else [str(point[0]), str(point[1]), str(point[2])]


def local_quotient(curve, points, prime: int, ell: int):
    field = GF(prime)
    reduced_curve = EllipticCurve(
        field,
        [field(a.numerator()) / field(a.denominator()) for a in curve.a_invariants()],
    )
    elements = list(reduced_curve)
    multiples = {ell * point for point in elements}
    unassigned = set(elements)
    cosets = []
    coset_of = {}
    while unassigned:
        representative = next(iter(unassigned))
        coset = {representative + h for h in multiples}
        index = len(cosets)
        for point in coset:
            coset_of[point] = index
        unassigned -= coset
        cosets.append(coset)

    zero = coset_of[reduced_curve(0)]
    representatives = [next(iter(coset)) for coset in cosets]
    span = {zero: ()}
    while len(span) < len(cosets):
        basis_index = next(i for i in range(len(cosets)) if i not in span)
        new_span = {}
        for index, vector in span.items():
            for coefficient in range(ell):
                target = representatives[index] + coefficient * representatives[basis_index]
                new_span[coset_of[target]] = vector + (coefficient,)
        span = new_span

    vectors = []
    for point in points:
        if point.is_zero():
            reduced = reduced_curve(0)
        else:
            x_coordinate, y_coordinate = point[0], point[1]
            if x_coordinate.denominator() % prime == 0 or y_coordinate.denominator() % prime == 0:
                raise ZeroDivisionError("point denominator")
            reduced = reduced_curve(
                field(x_coordinate.numerator()) / field(x_coordinate.denominator()),
                field(y_coordinate.numerator()) / field(y_coordinate.denominator()),
            )
        vectors.append(span[coset_of[reduced]])
    return vectors, int(reduced_curve.cardinality()), len(next(iter(span.values())))


def build_certificate(split_certificate: Path, prime_bound: int = 700) -> dict:
    split = json.loads(split_certificate.read_text())
    parameter = QQ(split["common_fibre_parameter"])
    roots = [QQ(value) for value in split["three_source_roots"]]

    ring = PolynomialRing(QQ, names=("X", "Y", "Z"))
    X, Y, Z = ring.gens()
    first = (
        3972 * X**3 + 8080 * X**2 * Y - 65622 * X**2 * Z
        + 31679 * X * Y**2 - 104467 * X * Y * Z - 232614 * X * Z**2
        + 24484 * Y**3 - 15556 * Y**2 * Z - 173688 * Y * Z**2
    )
    second = (
        33084 * X**3 + 44912 * X**2 * Y - 62778 * X**2 * Z
        + 24409 * X * Y**2 - 70613 * X * Y * Z - 138714 * X * Z**2
        - 36220 * Y**3 - 122924 * Y**2 * Z + 347376 * Y * Z**2
        + 1042128 * Z**3
    )
    cubic = first + parameter * second
    origin = [QQ(value) for value in ORIGIN]
    morphism = EllipticCurve_from_cubic(cubic, origin, morphism=True)
    curve = morphism.codomain()
    cubic_curve = morphism.domain()

    base = [morphism(cubic_curve([QQ(value) for value in point])) for point in BASE_POINTS[:8]]
    first_differences = []
    second_differences = []
    packet_records = []
    for root in roots:
        first_point, second_point, reference_point = (
            [QQ(value) for value in plane_point(root, kind)]
            for kind in ("C1", "C2", "Q")
        )
        for point in (first_point, second_point, reference_point):
            if cubic(*point) != 0:
                raise AssertionError(("point off cubic", root, point))
        mapped_first, mapped_second, mapped_reference = (
            morphism(cubic_curve(point))
            for point in (first_point, second_point, reference_point)
        )
        first_differences.append(mapped_first - mapped_reference)
        second_differences.append(mapped_second - mapped_reference)
        packet_records.append({
            "source_root": str(root),
            "C1": coords(mapped_first),
            "C2": coords(mapped_second),
            "Q": coords(mapped_reference),
            "D1": coords(first_differences[-1]),
            "D2": coords(second_differences[-1]),
        })

    all_points = base + [point for pair in zip(first_differences, second_differences) for point in pair]
    relation_1 = [-1, -1, -1, -1, -4, 5, 2, -1, 3, 0, 3, 0, 3, 0]
    relation_2 = [-4, 2, 2, -4, -1, 5, -1, 2, 0, 3, 0, 3, 0, 3]
    relation_checks = [
        sum((coefficient * point for coefficient, point in zip(relation, all_points)), curve(0)).is_zero()
        for relation in (relation_1, relation_2)
    ]
    if relation_checks != [True, True]:
        raise AssertionError(("trace relation failure", relation_checks))

    selected = base + [
        first_differences[0], second_differences[0],
        first_differences[1], second_differences[1],
    ]
    labels = [f"B{i}" for i in range(8)] + ["D1_0", "D2_0", "D1_1", "D2_1"]
    ell = 5
    rows = [[] for _ in selected]
    local_records = []
    skipped_primes = []
    torsion_witness = None
    rank = 0
    discriminant = curve.discriminant()
    for prime_element in prime_range(5, prime_bound):
        prime = int(prime_element)
        try:
            if discriminant.numerator() % prime == 0 or discriminant.denominator() % prime == 0:
                continue
            vectors, order, dimension = local_quotient(curve, selected, prime, ell)
        except Exception as exc:
            skipped_primes.append({"prime": prime, "error": repr(exc)})
            continue
        if torsion_witness is None and order % ell:
            torsion_witness = {"prime": prime, "group_order": order}
        if dimension:
            local_records.append({
                "prime": prime,
                "group_order": order,
                "quotient_dimension": dimension,
            })
            for index, vector in enumerate(vectors):
                rows[index].extend(vector)
            rank = int(matrix(GF(ell), rows).rank())
            if rank == len(selected) and torsion_witness is not None:
                break

    if rank != 12 or torsion_witness is None:
        raise AssertionError(("rank certificate failure", rank, torsion_witness))

    a_invariants = [str(value) for value in curve.a_invariants()]
    discriminant_text = str(discriminant)
    selected_coordinates = [coords(point) for point in selected]
    payload = {
        "schema_version": 1,
        "certificate_id": "small-split-e8-tangent-fibre-rank12-sage-v1",
        "claim_status": "independent SageMath rank-12 finite-reduction certificate",
        "sage_version": str(version()),
        "parent_split_certificate_sha256": split["record_sha256"],
        "fibre_parameter": str(parameter),
        "source_roots": [str(root) for root in roots],
        "weierstrass_a_invariants_sha256": canonical_hash(a_invariants),
        "weierstrass_discriminant_sha256": hashlib.sha256(discriminant_text.encode()).hexdigest(),
        "selected_labels": labels,
        "selected_points_sha256": canonical_hash(selected_coordinates),
        "packet_points_sha256": canonical_hash(packet_records),
        "trace_relation_checks": relation_checks,
        "ell": ell,
        "finite_reduction_rank": rank,
        "rows": rows,
        "local_records": local_records,
        "torsion_witness": torsion_witness,
        "proved_displayed_subgroup_rank_12": True,
        "truth_note": (
            "The finite-reduction matrix proves twelve independent points. "
            "Together with two exact independent trace relations among fourteen displayed points, "
            "the displayed subgroup has rank exactly twelve. This is not a rank-30 result."
        ),
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", type=Path)
    parser.add_argument("--prime-bound", type=int, default=700)
    args = parser.parse_args()

    payload = build_certificate(args.split_certificate, args.prime_bound)
    if args.check:
        previous = json.loads(args.check.read_text())
        if payload != previous:
            raise SystemExit("certificate mismatch")
        print("VERIFIED Sage rank-12 certificate sha256=" + payload["record_sha256"])
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "finite_reduction_rank": payload["finite_reduction_rank"],
        "used_primes": len(payload["local_records"]),
        "torsion_witness": payload["torsion_witness"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
