#!/usr/bin/env python3
"""Probe exact finite-field packet ranks with the official Magma calculator.

Magma's `AnalyticInformation` for elliptic curves over finite-field function
fields returns the arithmetic rank predicted by Artin--Tate; for rational and
K3 elliptic surfaces this rank is unconditional.  Each quadratic character is
submitted as an independent calculator request so a timeout cannot corrupt the
other channels.

The output remains finite-characteristic discovery evidence.  It is never
promoted to a characteristic-zero rank lower bound without explicit Q(t)
sections and an independent certificate.
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from fractions import Fraction
from pathlib import Path

SERVER = "http://magma.maths.usyd.edu.au/xml/calculator.xml"
REFERER = "http://magma.maths.usyd.edu.au/calc/"
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
        raise ZeroDivisionError(f"prime {p} divides coefficient denominator")
    return (value.numerator % p) * pow(value.denominator % p, -1, p) % p


def parse_support(token: str):
    token = token.strip().lower()
    if token in {"inf", "infinity", "oo"}:
        return None
    return int(token)


def label(value):
    return "infinity" if value is None else str(value)


def d_expression(a, b, p: int) -> str:
    if a is None:
        a, b = b, a
    aa = a % p
    if b is None:
        return f"t-F!{aa}"
    bb = b % p
    return f"(t-F!{aa})*(t-F!{bb})"


def magma_program(prime: int, a, b) -> str:
    c = {name: mod_fraction(value, prime) for name, value in RATIONALS.items()}
    d = d_expression(a, b, prime)
    return f'''SetColumns(0);
F := GF({prime});
K<t> := FunctionField(F);
a4 := K!(F!{c['p0']}) + K!(F!{c['p1']})*t + K!(F!{c['p2']})*t^2;
a6 := K!(F!{c['q0']}) + K!(F!{c['q1']})*t + K!(F!{c['q2']})*t^2 + K!(F!{c['q3']})*t^3 + K!(F!{c['q4']})*t^4 + t^5;
E := EllipticCurve([K|0,t^2,0,a4,a6]);
assert Discriminant(E) ne 0;
d := K!({d});
Et := QuadraticTwist(E,d);
assert Discriminant(Et) ne 0;
AI := AnalyticInformation(Et);
print "R30_MAGMA_VERSION", GetVersion();
print "R30_PRIME", {prime};
print "R30_CHARACTER", "{label(a)}|{label(b)}";
print "R30_ANALYTIC_INFORMATION", AI;
print "R30_ARITHMETIC_RANK", AI[1];
print "R30_GEOMETRIC_RANK", AI[2];
print "R30_PACKET_PROBE_PASS";
'''


def submit(name: str, code: str, output_dir: Path) -> dict:
    started = time.time()
    data = urllib.parse.urlencode({"input": code}).encode()
    request = urllib.request.Request(
        SERVER,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/xml,text/xml,text/html",
            "Referer": REFERER,
            "User-Agent": "elliptic-rank-30-packet-probe/1.0",
        },
        method="POST",
    )
    raw = b""
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read()
            http_status = response.status
        xml_path = output_dir / f"{name}.xml"
        txt_path = output_dir / f"{name}.txt"
        xml_path.write_bytes(raw)
        root = ET.fromstring(raw)
        warning = " ".join("".join(w.itertext()).strip() for w in root.findall(".//warning")).strip()
        lines = ["".join(line.itertext()) for line in root.findall(".//results/line")]
        output = "\n".join(lines)
        txt_path.write_text(output + "\n", encoding="utf-8")
        error_markers = (
            "Runtime error", "User error", "Internal error", "Assertion failed",
            "Identifier '", "Syntax error", "time limit",
        )
        has_error = any(marker.lower() in output.lower() for marker in error_markers)
        passed = "R30_PACKET_PROBE_PASS" in output and not warning and not has_error
        rank = None
        geometric_rank = None
        for line in output.splitlines():
            fields = line.strip().split()
            if fields[:1] == ["R30_ARITHMETIC_RANK"] and len(fields) >= 2:
                rank = int(fields[1])
            if fields[:1] == ["R30_GEOMETRIC_RANK"] and len(fields) >= 2:
                geometric_rank = int(fields[1])
        return {
            "name": name,
            "status": "pass" if passed else "fail",
            "rank": rank,
            "geometric_rank": geometric_rank,
            "warning": warning,
            "has_error_marker": has_error,
            "http_status": http_status,
            "elapsed_seconds": time.time() - started,
            "output_file": txt_path.name,
            "xml_file": xml_path.name,
            "output_tail": output[-3000:],
        }
    except Exception:
        if raw:
            (output_dir / f"{name}.xml").write_bytes(raw)
        return {
            "name": name,
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--triple", default="0,1,inf")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    support = [parse_support(token) for token in args.triple.split(",")]
    if len(support) != 3 or len({label(x) for x in support}) != 3:
        raise ValueError("--triple must contain three distinct branch points")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for a, b in ((support[0], support[1]), (support[0], support[2]), (support[1], support[2])):
        name = f"magma_packet_p{args.prime}_{label(a)}_{label(b)}".replace("-", "m")
        code = magma_program(args.prime, a, b)
        (args.output_dir / f"{name}.m").write_text(code, encoding="utf-8")
        results.append(submit(name, code, args.output_dir))

    ranks = [row["rank"] for row in results if row["status"] == "pass" and row["rank"] is not None]
    summary = {
        "status": "pass" if len(ranks) == 3 else "incomplete",
        "prime": args.prime,
        "support": [label(x) for x in support],
        "results": results,
        "twist_rank_sum": sum(ranks) if len(ranks) == 3 else None,
        "packet_score_using_characteristic_zero_base_rank_8": 8 + sum(ranks) if len(ranks) == 3 else None,
        "truth_note": (
            "exact finite-characteristic rank data only; the corresponding "
            "characteristic-zero twists require explicit Q(t)-sections"
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"],
        "prime": args.prime,
        "support": summary["support"],
        "ranks": ranks,
        "packet_score": summary["packet_score_using_characteristic_zero_base_rank_8"],
    }, sort_keys=True))
    # Preserve timeout/failure evidence; the JSON status is authoritative.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
