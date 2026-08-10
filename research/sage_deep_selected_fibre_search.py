from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import signal
import traceback
from pathlib import Path

from sage.all import EllipticCurve, GF, QQ, matrix, prime_range, version


class SearchTimeout(Exception):
    pass


def alarm_handler(signum, frame):
    raise SearchTimeout("bounded search timed out")


def with_timeout(seconds, function, *args, **kwargs):
    previous = signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(seconds)
    try:
        return function(*args, **kwargs), None
    except Exception as exc:
        return None, repr(exc)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def parse_point(E, values):
    if values == ["0"] or values == [0]:
        return E(0)
    if len(values) == 2:
        return E(QQ(values[0]), QQ(values[1]))
    if len(values) == 3:
        return E([QQ(values[0]), QQ(values[1]), QQ(values[2])])
    raise ValueError(values)


def coordinates(point):
    return ["0"] if point.is_zero() else [str(point[0]), str(point[1]), str(point[2])]


def rank_mod(rows, ell):
    return 0 if not rows else int(matrix(GF(ell), rows).rank())


def local_quotient(E, points, prime, ell):
    field = GF(prime)
    reduced = EllipticCurve(
        field,
        [field(a.numerator()) / field(a.denominator()) for a in E.a_invariants()],
    )
    elements = list(reduced)
    multiples = {ell * point for point in elements}
    unassigned = set(elements)
    cosets = []
    which = {}
    while unassigned:
        representative = next(iter(unassigned))
        coset = {representative + multiple for multiple in multiples}
        index = len(cosets)
        for point in coset:
            which[point] = index
        unassigned -= coset
        cosets.append(coset)

    zero = which[reduced(0)]
    representatives = [next(iter(coset)) for coset in cosets]
    coordinates_by_coset = {zero: ()}
    while len(coordinates_by_coset) < len(cosets):
        basis = next(i for i in range(len(cosets)) if i not in coordinates_by_coset)
        new = {}
        for index, vector in coordinates_by_coset.items():
            for coefficient in range(ell):
                target = which[representatives[index] + coefficient * representatives[basis]]
                new[target] = vector + (coefficient,)
        coordinates_by_coset = new

    vectors = []
    for point in points:
        if point.is_zero():
            reduced_point = reduced(0)
        else:
            x_coordinate, y_coordinate = point[0], point[1]
            if x_coordinate.denominator() % prime == 0 or y_coordinate.denominator() % prime == 0:
                raise ZeroDivisionError
            reduced_point = reduced(
                field(x_coordinate.numerator()) / field(x_coordinate.denominator()),
                field(y_coordinate.numerator()) / field(y_coordinate.denominator()),
            )
        vectors.append(coordinates_by_coset[which[reduced_point]])
    return vectors, int(reduced.cardinality()), len(next(iter(coordinates_by_coset.values())))


def finite_reduction_certificate(E, points, prime_bound):
    discriminant = E.discriminant()
    for ell in (2, 3, 5, 7, 11, 13, 17):
        rows = [[] for _ in points]
        local_records = []
        torsion_witness = None
        for prime0 in prime_range(5, prime_bound):
            prime = int(prime0)
            try:
                if discriminant.numerator() % prime == 0 or discriminant.denominator() % prime == 0:
                    continue
                vectors, order, dimension = local_quotient(E, points, prime, ell)
            except Exception:
                continue
            if torsion_witness is None and order % ell:
                torsion_witness = {"prime": prime, "group_order": order}
            if dimension:
                for index, vector in enumerate(vectors):
                    rows[index].extend(vector)
                local_records.append(
                    {"prime": prime, "group_order": order, "quotient_dimension": dimension}
                )
                rank = rank_mod(rows, ell)
                if rank == len(points) and torsion_witness is not None:
                    return {
                        "ell": ell,
                        "rank": rank,
                        "rows": rows,
                        "local_records": local_records,
                        "torsion_witness": torsion_witness,
                    }
    return None


def select_records(payload, count):
    records = list(payload.get("records") or [])
    records.sort(
        key=lambda record: (
            0 if "timed out" in str(record.get("sage_generator_search_error") or "").lower() else 1,
            -(int(record["sage_rank_estimate_unproved"]) if record.get("sage_rank_estimate_unproved") is not None else -1),
            int(record.get("coefficient_height_bits") or 10**9),
        )
    )
    return records[:count]


def append_points(target, candidates):
    for candidate in candidates or []:
        try:
            point = candidate
            if isinstance(candidate, (list, tuple)) and candidate and not hasattr(candidate, "curve"):
                point = target[0].curve()(candidate)
            if point.is_zero():
                continue
            if any(point == known or point == -known for known in target):
                continue
            target.append(point)
        except Exception:
            continue


def search_curve(record, timeout_seconds, prime_bound):
    E = EllipticCurve(QQ, [QQ(value) for value in record["a_invariants"]])
    known = [parse_point(E, values) for values in record["known_rank12_points"]]
    baseline = finite_reduction_certificate(E, known, prime_bound)
    if baseline is None or baseline["rank"] != 12:
        raise AssertionError("the displayed rank-12 subgroup did not replay")

    pool = []
    attempts = []

    methods = []
    methods.append(("gens", lambda: E.gens(proof=False)))
    methods.append(("gens_rank13", lambda: E.gens(proof=False, rank=13)))
    if hasattr(E, "point_search"):
        methods.extend(
            [
                ("point_search_18", lambda: E.point_search(18)),
                ("point_search_22", lambda: E.point_search(22)),
            ]
        )
    if hasattr(E, "integral_points"):
        methods.append(("integral_points", lambda: E.integral_points(both_signs=True)))

    for name, method in methods:
        value, error = with_timeout(timeout_seconds, method)
        attempts.append(
            {
                "method": name,
                "error": error,
                "result_count": None if value is None else len(value),
            }
        )
        if value is not None:
            append_points(pool, value)

    current = list(known)
    accepted = []
    for point in pool:
        if any(point == known_point or point == -known_point for known_point in current):
            continue
        certificate = finite_reduction_certificate(E, current + [point], prime_bound)
        if certificate is not None and certificate["rank"] == len(current) + 1:
            accepted.append({"point": coordinates(point), "certificate": certificate})
            current.append(point)

    return {
        "multiple": record.get("multiple"),
        "fibre_parameter": record.get("fibre_parameter"),
        "a_invariants": record.get("a_invariants"),
        "coefficient_height_bits": record.get("coefficient_height_bits"),
        "baseline_certificate": baseline,
        "attempts": attempts,
        "candidate_pool": [coordinates(point) for point in pool],
        "accepted_new_generators": accepted,
        "certified_displayed_rank": len(current),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--curves", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--prime-bound", type=int, default=1000)
    args = parser.parse_args()

    try:
        source = json.loads(args.input.read_text())
        selected = select_records(source, args.curves)
        records = []
        for record in selected:
            try:
                records.append(search_curve(record, args.timeout, args.prime_bound))
            except Exception as exc:
                records.append(
                    {
                        "multiple": record.get("multiple"),
                        "status": "error",
                        "error": repr(exc),
                        "traceback": traceback.format_exc(),
                        "certified_displayed_rank": 0,
                    }
                )
        payload = {
            "schema_version": 1,
            "status": "completed",
            "sage_version": str(version()),
            "source_record_sha256": source.get("record_sha256"),
            "curves_searched": len(records),
            "records": records,
            "maximum_certified_displayed_rank": max(
                (record.get("certified_displayed_rank", 0) for record in records), default=0
            ),
            "winner_count": sum(
                bool(record.get("accepted_new_generators")) for record in records
            ),
            "solved_rank30": any(
                record.get("certified_displayed_rank", 0) >= 30 for record in records
            ),
            "truth_note": (
                "Only generators accompanied by a full finite-reduction independence certificate "
                "increase the certified displayed rank."
            ),
        }
    except Exception as exc:
        payload = {
            "schema_version": 1,
            "status": "error",
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "solved_rank30": False,
        }

    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": payload.get("status"),
                "curves_searched": payload.get("curves_searched"),
                "winner_count": payload.get("winner_count"),
                "maximum_certified_displayed_rank": payload.get(
                    "maximum_certified_displayed_rank"
                ),
                "solved_rank30": payload.get("solved_rank30"),
                "error": payload.get("error"),
                "record_sha256": payload.get("record_sha256"),
            },
            indent=2,
            sort_keys=True,
        )
    )
    if payload.get("status") == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
