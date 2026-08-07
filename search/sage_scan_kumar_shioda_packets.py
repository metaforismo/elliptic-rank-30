#!/usr/bin/env sage-python
"""Finite-field rank scan for V4 packets on the certified Kumar--Shioda E8 surface.

The characteristic-zero base surface has exact rank 8 over Q(t). For a fixed
good prime, this script reduces the surface modulo p and computes Sage's
function-field ranks for every quadratic twist supported on a pair from a
small set of rational branch points. It then scores every three-branch V4
packet by

    8 + rank(E^d_ab) + rank(E^d_ac) + rank(E^d_bc).

The output is a discovery and upper-sieve laboratory only. A high rank after
reduction is not a characteristic-zero lower bound. Promotion requires
explicit Q(t)-sections and an independent certificate.
"""
from __future__ import annotations

import argparse
import itertools
import json
import traceback
from fractions import Fraction
from pathlib import Path

from sage.all import GF, FunctionField, EllipticCurve

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


def reduce_fraction(value: Fraction, field):
    den = field(value.denominator)
    if den == 0:
        raise ZeroDivisionError(f"prime divides denominator {value.denominator}")
    return field(value.numerator) / den


def support_label(value):
    return "infinity" if value is None else str(value)


def pair_key(a, b):
    return tuple(sorted((support_label(a), support_label(b))))


def twist_parameter(t, field, a, b):
    if a is None:
        a, b = b, a
    aa = field(a)
    if b is None:
        return t - aa
    return (t - aa) * (t - field(b))


def build_surface(prime):
    field = GF(prime)
    K = FunctionField(field, "t")
    t = K.gen()
    c = {name: reduce_fraction(value, field) for name, value in RATIONALS.items()}
    a4 = c["p0"] + c["p1"] * t + c["p2"] * t**2
    a6 = c["q0"] + c["q1"] * t + c["q2"] * t**2 + c["q3"] * t**3 + c["q4"] * t**4 + t**5
    E = EllipticCurve(K, [0, t**2, 0, a4, a6])
    if E.discriminant() == 0:
        raise RuntimeError("singular reduction")
    return field, K, t, E


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--bound", type=int, default=4,
                        help="use finite branch points -bound,...,bound plus infinity")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = {
        "status": "started",
        "prime": args.prime,
        "bound": args.bound,
        "truth_note": (
            "finite-characteristic discovery/upper-sieve evidence only; "
            "no Q(t) or rank-30 lower-bound claim"
        ),
    }

    try:
        field, K, t, E = build_surface(args.prime)
        output["reduction_base_rank"] = int(E.rank())
        output["reduction_discriminant"] = str(E.discriminant())
        supports = list(range(-args.bound, args.bound + 1)) + [None]
        if len({field(x) for x in supports if x is not None}) != len(supports) - 1:
            raise ValueError("branch support collides modulo the chosen prime")
        output["supports"] = [support_label(x) for x in supports]

        pair_ranks = {}
        for a, b in itertools.combinations(supports, 2):
            key = pair_key(a, b)
            d = twist_parameter(t, field, a, b)
            record = {"d": str(d)}
            try:
                twist = E.quadratic_twist(d)
                record["rank"] = int(twist.rank())
                record["discriminant"] = str(twist.discriminant())
            except Exception as exc:
                record["error"] = repr(exc)
                record["traceback"] = traceback.format_exc()
            pair_ranks["|".join(key)] = record
        output["pair_ranks"] = pair_ranks

        packets = []
        for triple in itertools.combinations(supports, 3):
            ranks = []
            keys = []
            ok = True
            for a, b in itertools.combinations(triple, 2):
                key = pair_key(a, b)
                rec = pair_ranks["|".join(key)]
                keys.append("|".join(key))
                if "rank" not in rec:
                    ok = False
                    break
                ranks.append(rec["rank"])
            if ok:
                packets.append({
                    "support": [support_label(x) for x in triple],
                    "pair_keys": keys,
                    "twist_ranks": ranks,
                    "characteristic_zero_target_score": 8 + sum(ranks),
                    "reduction_total_score": output["reduction_base_rank"] + sum(ranks),
                })
        packets.sort(
            key=lambda row: (
                row["characteristic_zero_target_score"],
                sorted(row["twist_ranks"], reverse=True),
            ),
            reverse=True,
        )
        output["packet_count"] = len(packets)
        output["top_packets"] = packets[:25]
        output["best_target_score"] = packets[0]["characteristic_zero_target_score"] if packets else None
        output["status"] = "pass"
    except Exception as exc:
        output["status"] = "error"
        output["error"] = repr(exc)
        output["traceback"] = traceback.format_exc()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": output["status"],
        "prime": args.prime,
        "reduction_base_rank": output.get("reduction_base_rank"),
        "best_target_score": output.get("best_target_score"),
        "packet_count": output.get("packet_count"),
    }, sort_keys=True))
    return 0 if output["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
