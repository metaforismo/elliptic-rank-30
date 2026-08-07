#!/usr/bin/env python3
"""Scan all nondegenerate markings modulo p for maximality of the k20 K3 channel.

The channel y^2=x^3+v^4*P_mu(v) has geometric rank capacity four.  A rational
marking capable of characteristic-zero rank 30 must reduce to geometric rank
four at every good prime.  This script computes the exact finite-field
geometric rank with Magma's AnalyticInformation and returns the surviving
residue classes.  It is a necessary-condition sieve only.
"""
from __future__ import annotations

import argparse
import json
import time
import traceback
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

SERVER = "http://magma.maths.usyd.edu.au/xml/calculator.xml"
REFERER = "http://magma.maths.usyd.edu.au/calc/"


def program(prime: int, mu: int) -> str:
    return f'''SetColumns(0);
F:=GF({prime}); K<v>:=FunctionField(F);
mu:=F!({mu});
a:=mu^2-mu+1; rmu:=mu*(mu-1);
assert rmu ne 0;
b:=a^3/rmu;
P:=K!(b^2)*v^2-K!(a^3)*(v+1)^3;
E:=EllipticCurve([K|0,0,0,0,v^4*P]);
assert Discriminant(E) ne 0;
AI:=AnalyticInformation(E);
print "R30_PRIME",{prime};
print "R30_MU",{mu};
print "R30_ARITHMETIC_RANK",AI[1];
print "R30_GEOMETRIC_RANK",AI[2];
print "R30_K20_PASS";
'''


def submit(prime: int, mu: int, outdir: Path) -> dict:
    name = f"k20_p{prime}_mu{mu}"
    code = program(prime, mu)
    (outdir / f"{name}.m").write_text(code)
    started = time.time()
    data = urllib.parse.urlencode({"input": code}).encode()
    request = urllib.request.Request(
        SERVER,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/xml,text/xml,text/html",
            "Referer": REFERER,
            "User-Agent": "elliptic-rank-30-k20-residue-scan/1.0",
        },
        method="POST",
    )
    raw = b""
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read()
            http = response.status
        (outdir / f"{name}.xml").write_bytes(raw)
        root = ET.fromstring(raw)
        warning = " ".join(
            "".join(node.itertext()).strip() for node in root.findall(".//warning")
        ).strip()
        lines = ["".join(node.itertext()) for node in root.findall(".//results/line")]
        text = "\n".join(lines)
        (outdir / f"{name}.txt").write_text(text + "\n")
        arithmetic = geometric = None
        for line in text.splitlines():
            fields = line.split()
            if fields[:1] == ["R30_ARITHMETIC_RANK"] and len(fields) > 1:
                arithmetic = int(fields[1])
            if fields[:1] == ["R30_GEOMETRIC_RANK"] and len(fields) > 1:
                geometric = int(fields[1])
        passed = "R30_K20_PASS" in text and not warning
        return {
            "mu": mu,
            "status": "pass" if passed else "fail",
            "arithmetic_rank": arithmetic,
            "geometric_rank": geometric,
            "warning": warning,
            "http_status": http,
            "elapsed_seconds": time.time() - started,
        }
    except Exception:
        if raw:
            (outdir / f"{name}.xml").write_bytes(raw)
        return {
            "mu": mu,
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    p = args.prime
    args.output_dir.mkdir(parents=True, exist_ok=True)
    excluded = {0, 1, 2, p - 1, pow(2, p - 2, p)}
    residues = [mu for mu in range(p) if mu not in excluded]
    rows = [submit(p, mu, args.output_dir) for mu in residues]
    maximal = [row["mu"] for row in rows if row["status"] == "pass" and row["geometric_rank"] == 4]
    rejected = [row["mu"] for row in rows if row["status"] == "pass" and row["geometric_rank"] < 4]
    result = {
        "status": "pass",
        "prime": p,
        "excluded_degenerate_residues": sorted(excluded),
        "rows": rows,
        "maximal_residues": maximal,
        "rejected_residues": rejected,
        "incomplete_residues": [row["mu"] for row in rows if row["status"] != "pass"],
        "truth_status": "exact finite-characteristic necessary-condition sieve",
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "prime": p,
        "maximal_residues": maximal,
        "rejected_count": len(rejected),
        "incomplete_count": len(result["incomplete_residues"]),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
