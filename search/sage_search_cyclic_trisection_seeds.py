#!/usr/bin/env sage-python
"""Search exact finite-field cyclic trisection seeds on the split E8 surface.

The search uses a height-six trace section S=alpha+beta and the generic-fibre
basis

    1, x, g_S=(y+y_S)/(x-x_S)

of L(2O+S). For

    f=a(t)+b(t)x+g_S,

the x-coordinates of its three zeros satisfy an explicit cubic. An
irreducible cubic with square discriminant defines a cyclic cubic extension.
Finite-characteristic hits are lifting seeds only.
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import traceback
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    FunctionField,
    PolynomialRing,
    QQ,
    Matrix,
    sage_eval,
)

ROOT = Path(__file__).resolve().parents[1]
NORMALIZER = ROOT / "search" / "normalize_kumar_shioda_certificate.py"
SOURCE = ROOT / "certificates" / "kumar_shioda_split_e8_certificate.json"


def parse_expression(text, t):
    return sage_eval(str(text).replace("^", "**"), locals={"t": t})


def ensure_normalized(path):
    if path.exists():
        return
    subprocess.run(
        [
            sys.executable,
            str(NORMALIZER),
            "--source",
            str(SOURCE),
            "--output",
            str(path),
        ],
        check=True,
        cwd=ROOT,
    )


def choose_root_pair(height_matrix, points):
    for i in range(8):
        if height_matrix[i, i] != 2:
            continue
        for j in range(i + 1, 8):
            if height_matrix[j, j] != 2:
                continue
            pairing = height_matrix[i, j]
            if pairing == 1:
                return {
                    "indices": [i, j],
                    "signs": [1, 1],
                    "alpha": points[i],
                    "beta": points[j],
                }
            if pairing == -1:
                return {
                    "indices": [i, j],
                    "signs": [1, -1],
                    "alpha": points[i],
                    "beta": -points[j],
                }
    raise RuntimeError("the certified basis contains no immediately visible root pair")


def rational_function_is_square(value):
    try:
        return bool(value.is_square())
    except Exception:
        numerator = value.numerator()
        denominator = value.denominator()
        try:
            return bool(numerator.is_square() and denominator.is_square())
        except Exception:
            return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--max-hits", type=int, default=20)
    parser.add_argument("--normalized", type=Path, default=ROOT / "search" / "kumar_shioda_normalized.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = {
        "status": "started",
        "prime": args.prime,
        "truth_status": "finite-characteristic lifting seeds only; no Q(t) or rank-30 claim",
    }
    try:
        ensure_normalized(args.normalized)
        data = json.loads(args.normalized.read_text())

        field = GF(args.prime)
        function_field = FunctionField(field, "t")
        t = function_field.gen()
        a4 = function_field(parse_expression(data["weierstrass_model"]["a4"], t))
        a6 = function_field(parse_expression(data["weierstrass_model"]["a6"], t))
        curve = EllipticCurve(function_field, [0, 0, 0, a4, a6])

        points = []
        for record in data["sections"]:
            x = function_field(parse_expression(record["x"], t))
            y = function_field(parse_expression(record["y"], t))
            point = curve([x, y])
            points.append(point)
        assert len(points) == 8

        height_matrix = Matrix(QQ, [[QQ(entry) for entry in row] for row in data["height_matrix"]])
        pair = choose_root_pair(height_matrix, points)
        trace = pair["alpha"] + pair["beta"]
        if trace.is_zero():
            raise RuntimeError("chosen height-six trace section is zero")
        u = trace[0]
        v = trace[1]

        polynomial_ring = PolynomialRing(function_field, "X")
        X = polynomial_ring.gen()
        hits = []
        checked = 0
        for a0, a1, b0, b1 in itertools.product(field, repeat=4):
            if b0 == 0 and b1 == 0:
                continue
            checked += 1
            a = function_field(a0 + a1 * t)
            b = function_field(b0 + b1 * t)
            cubic = polynomial_ring(
                (a + b * X) ** 2 * (X - u)
                + 2 * v * (a + b * X)
                - (X ** 2 + u * X + u ** 2 + a4)
            )
            if cubic.degree() != 3 or cubic.discriminant() == 0:
                continue
            if not cubic.is_irreducible():
                continue
            discriminant = cubic.discriminant()
            if not rational_function_is_square(discriminant):
                continue

            hit = {
                "a0": int(a0),
                "a1": int(a1),
                "b0": int(b0),
                "b1": int(b1),
                "a": str(a),
                "b": str(b),
                "trace_x": str(u),
                "trace_y": str(v),
                "cubic": str(cubic),
                "cubic_discriminant": str(discriminant),
                "irreducible": True,
                "square_discriminant": True,
            }
            try:
                extension = function_field.extension(cubic, names=("z",))
                hit["extension_genus"] = int(extension.genus())
            except Exception as exc:
                hit["genus_error"] = repr(exc)
            hits.append(hit)
            if len(hits) >= args.max_hits:
                break

        result.update(
            {
                "status": "pass",
                "checked_parameter_tuples": checked,
                "root_pair_indices": pair["indices"],
                "root_pair_signs": pair["signs"],
                "hit_count": len(hits),
                "hits": hits,
            }
        )
    except Exception as exc:
        result["status"] = "error"
        result["error"] = repr(exc)
        result["traceback"] = traceback.format_exc()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "prime": args.prime,
                "checked": result.get("checked_parameter_tuples"),
                "hits": result.get("hit_count"),
                "error": result.get("error"),
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
