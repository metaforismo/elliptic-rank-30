#!/usr/bin/env python3
"""Finite-field Picard/rank sieve for the three K3 channels of the marked j=0 packet.

For a rational marking mu, put a=mu^2-mu+1, b=a^3/(mu(mu-1)),
P(v)=b^2*v^2-a^3*(v+1)^3 and D(v)=1+4*v.  The full C2 x C3 packet has six
channels.  Three are rational elliptic surfaces and attain their geometric
capacity automatically.  Rank 30 therefore requires the three K3 channels

  k20: v^4*P,
  k11: v^2*D^3*P,
  k21: v^4*D^3*P

to have geometric ranks 4, 6 and 6.  Magma's AnalyticInformation over finite
function fields gives an exact discovery/rejection sieve.  No finite-field
result is promoted to characteristic zero without explicit sections.
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
CHANNELS = {
    "k20": (4, 0, 4),
    "k11": (2, 3, 6),
    "k21": (4, 3, 6),
}


def magma_fraction(token: str) -> str:
    value = Fraction(token.strip())
    return f"F!({value.numerator})/F!({value.denominator})"


def program(prime: int, mu_token: str, name: str) -> str:
    v_power, d_power, target = CHANNELS[name]
    mu_expr = magma_fraction(mu_token)
    return f'''SetColumns(0);
F:=GF({prime}); K<v>:=FunctionField(F);
mu:={mu_expr};
assert mu ne 0 and mu ne 1;
a:=mu^2-mu+1; rmu:=mu*(mu-1); b:=a^3/rmu;
P:=K!(b^2)*v^2-K!(a^3)*(v+1)^3;
D:=1+4*v;
A6:=v^{v_power}*D^{d_power}*P;
assert A6 ne 0;
E:=EllipticCurve([K|0,0,0,0,A6]);
assert Discriminant(E) ne 0;
AI:=AnalyticInformation(E);
print "R30_PRIME",{prime};
print "R30_MU",mu;
print "R30_CHANNEL","{name}";
print "R30_TARGET",{target};
print "R30_ARITHMETIC_RANK",AI[1];
print "R30_GEOMETRIC_RANK",AI[2];
print "R30_SIX_CHANNEL_PASS";
'''


def submit(name: str, code: str, outdir: Path) -> dict:
    started = time.time()
    data = urllib.parse.urlencode({"input": code}).encode()
    request = urllib.request.Request(
        SERVER,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/xml,text/xml,text/html",
            "Referer": REFERER,
            "User-Agent": "elliptic-rank-30-six-channel/1.0",
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
        passed = "R30_SIX_CHANNEL_PASS" in text and not warning
        return {
            "status": "pass" if passed else "fail",
            "arithmetic_rank": arithmetic,
            "geometric_rank": geometric,
            "warning": warning,
            "http_status": http,
            "elapsed_seconds": time.time() - started,
            "output_tail": text[-2500:],
        }
    except Exception:
        if raw:
            (outdir / f"{name}.xml").write_bytes(raw)
        return {
            "status": "exception",
            "elapsed_seconds": time.time() - started,
            "traceback": traceback.format_exc(),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--mu", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for channel, (_, _, target) in CHANNELS.items():
        name = f"six_channel_p{args.prime}_mu{args.mu.replace('/', '_')}_{channel}"
        code = program(args.prime, args.mu, channel)
        (args.output_dir / f"{name}.m").write_text(code)
        row = submit(name, code, args.output_dir)
        row.update({
            "prime": args.prime,
            "mu": args.mu,
            "channel": channel,
            "target_geometric_rank": target,
        })
        rows.append(row)

    complete = all(row["status"] == "pass" for row in rows)
    maximal = complete and all(
        row["geometric_rank"] == row["target_geometric_rank"] for row in rows
    )
    rejected = any(
        row["status"] == "pass"
        and row["geometric_rank"] < row["target_geometric_rank"]
        for row in rows
    )
    result = {
        "status": "pass",
        "prime": args.prime,
        "mu": args.mu,
        "channels": rows,
        "all_channels_completed": complete,
        "all_three_k3_channels_maximal": maximal,
        "characteristic_zero_candidate_rejected_by_this_prime": rejected,
        "truth_status": (
            "exact finite-characteristic rank sieve; characteristic-zero rank "
            "requires explicit sections and independent certification"
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "prime": args.prime,
        "mu": args.mu,
        "ranks": [row.get("geometric_rank") for row in rows],
        "maximal": maximal,
        "rejected": rejected,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
