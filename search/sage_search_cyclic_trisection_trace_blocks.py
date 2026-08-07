#!/usr/bin/env sage-python
"""Block search over all 6720 minimal E8 trace classes.

For each selected norm-six trace vector, construct the corresponding section
on the certified split-E8 rational elliptic surface and search exact finite-
field functions in L(2O+S). Irreducible square-discriminant cubics of genus
one, with no detected branch collision against the elliptic discriminant, are
retained as characteristic-zero lifting seeds only.
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

from sage.all import EllipticCurve, GF, FunctionField, Matrix, PolynomialRing, QQ, sage_eval

ROOT = Path(__file__).resolve().parents[1]
NORMALIZER = ROOT / "search" / "normalize_kumar_shioda_certificate.py"
SOURCE = ROOT / "certificates" / "kumar_shioda_split_e8_certificate.json"


def parse_expression(text, t):
    return sage_eval(str(text).replace("^", "**"), locals={"t": t})


def ensure_normalized(path):
    if path.exists():
        return
    subprocess.run(
        [sys.executable, str(NORMALIZER), "--source", str(SOURCE), "--output", str(path)],
        cwd=ROOT,
        check=True,
    )


def pairing(left, right, gram):
    return sum(
        left[i] * gram[i, j] * right[j]
        for i in range(8)
        for j in range(8)
    )


def enumerate_e8_roots(gram):
    simple = []
    for i in range(8):
        vector = tuple(1 if i == j else 0 for j in range(8))
        if pairing(vector, vector, gram) == 2:
            simple.append(vector)
    if len(simple) < 8:
        raise RuntimeError("the normalized basis is not an E8 root basis")

    roots = set(simple)
    roots.update(tuple(-entry for entry in root) for root in simple)
    queue = list(roots)
    while queue:
        vector = queue.pop()
        for root in simple:
            coefficient = pairing(vector, root, gram)
            reflected = tuple(vector[i] - coefficient * root[i] for i in range(8))
            if reflected not in roots:
                if pairing(reflected, reflected, gram) != 2:
                    raise AssertionError("Weyl reflection changed the root norm")
                roots.add(reflected)
                queue.append(reflected)
    if len(roots) != 240:
        raise RuntimeError(f"expected 240 E8 roots, obtained {len(roots)}")
    return sorted(roots)


def enumerate_norm_six_traces(roots, gram):
    traces = {}
    for i, alpha in enumerate(roots):
        for beta in roots[i + 1 :]:
            if pairing(alpha, beta, gram) != 1:
                continue
            trace = tuple(alpha[j] + beta[j] for j in range(8))
            if pairing(trace, trace, gram) != 6:
                raise AssertionError("root pair did not produce norm six")
            pair = (alpha, beta)
            if trace in traces and traces[trace] != pair:
                raise AssertionError("norm-six trace has more than one root-pair decomposition")
            traces[trace] = pair
    if len(traces) != 6720:
        raise RuntimeError(f"expected 6720 norm-six traces, obtained {len(traces)}")
    return sorted(traces), traces


def combine_point(curve, coefficients, basis_points):
    point = curve(0)
    for coefficient, basis_point in zip(coefficients, basis_points):
        if coefficient:
            point += int(coefficient) * basis_point
    return point


def is_square_rational_function(value):
    try:
        return bool(value.is_square())
    except Exception:
        return bool(value.numerator().is_square() and value.denominator().is_square())


def coefficient_rows(field, t, pattern):
    elements = list(field)
    if pattern == "linear-linear":
        for a0, a1, b0, b1 in itertools.product(elements, repeat=4):
            if b0 == 0 and b1 == 0:
                continue
            yield (a0, a1, b0, b1), a0 + a1 * t, b0 + b1 * t
    elif pattern == "quadratic-constant":
        for a0, a1, a2, b0 in itertools.product(elements, repeat=4):
            if b0 == 0:
                continue
            yield (a0, a1, a2, b0), a0 + a1 * t + a2 * t**2, b0
    elif pattern == "constant-quadratic":
        for a0, b0, b1, b2 in itertools.product(elements, repeat=4):
            if b0 == 0 and b1 == 0 and b2 == 0:
                continue
            yield (a0, b0, b1, b2), a0, b0 + b1 * t + b2 * t**2
    else:
        raise ValueError(f"unknown coefficient pattern: {pattern}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--trace-start", type=int, default=0)
    parser.add_argument("--trace-count", type=int, default=20)
    parser.add_argument(
        "--pattern",
        choices=("linear-linear", "quadratic-constant", "constant-quadratic"),
        default="linear-linear",
    )
    parser.add_argument("--max-hits", type=int, default=25)
    parser.add_argument(
        "--normalized",
        type=Path,
        default=ROOT / "search" / "kumar_shioda_normalized.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    started = time.time()
    output = {
        "status": "started",
        "prime": args.prime,
        "trace_start": args.trace_start,
        "trace_count": args.trace_count,
        "pattern": args.pattern,
        "truth_status": "finite-characteristic lifting seeds only; no Q(t) or rank-30 claim",
    }
    try:
        ensure_normalized(args.normalized)
        data = json.loads(args.normalized.read_text())
        gram = Matrix(QQ, [[QQ(entry) for entry in row] for row in data["height_matrix"]])
        roots = enumerate_e8_roots(gram)
        trace_vectors, root_pairs = enumerate_norm_six_traces(roots, gram)

        field = GF(args.prime)
        function_field = FunctionField(field, "t")
        t = function_field.gen()
        a4 = function_field(parse_expression(data["weierstrass_model"]["a4"], t))
        a6 = function_field(parse_expression(data["weierstrass_model"]["a6"], t))
        curve = EllipticCurve(function_field, [0, 0, 0, a4, a6])
        elliptic_delta = curve.discriminant()
        elliptic_bad = elliptic_delta.numerator() * elliptic_delta.denominator()

        basis_points = []
        for record in data["sections"]:
            x = function_field(parse_expression(record["x"], t))
            y = function_field(parse_expression(record["y"], t))
            basis_points.append(curve([x, y]))

        polynomial_ring = PolynomialRing(function_field, "X")
        X = polynomial_ring.gen()
        selected = trace_vectors[args.trace_start : args.trace_start + args.trace_count]
        hits = []
        trace_reports = []
        for local_index, trace_coefficients in enumerate(selected):
            global_index = args.trace_start + local_index
            trace_point = combine_point(curve, trace_coefficients, basis_points)
            if trace_point.is_zero():
                raise AssertionError("a norm-six trace coefficient vector mapped to zero")
            u = trace_point[0]
            v = trace_point[1]
            report = {
                "trace_index": global_index,
                "trace_coefficients": list(trace_coefficients),
                "root_pair_coefficients": [
                    list(root_pairs[trace_coefficients][0]),
                    list(root_pairs[trace_coefficients][1]),
                ],
                "trace_x": str(u),
                "trace_y": str(v),
                "parameter_tuples_checked": 0,
                "hit_count": 0,
            }
            for parameters, a_value, b_value in coefficient_rows(field, t, args.pattern):
                report["parameter_tuples_checked"] += 1
                a_value = function_field(a_value)
                b_value = function_field(b_value)
                cubic = polynomial_ring(
                    (a_value + b_value * X) ** 2 * (X - u)
                    + 2 * v * (a_value + b_value * X)
                    - (X**2 + u * X + u**2 + a4)
                )
                if cubic.degree() != 3:
                    continue
                discriminant = cubic.discriminant()
                if discriminant == 0 or not is_square_rational_function(discriminant):
                    continue
                if not cubic.is_irreducible():
                    continue

                hit = {
                    "trace_index": global_index,
                    "trace_coefficients": list(trace_coefficients),
                    "parameters": [int(value) for value in parameters],
                    "a": str(a_value),
                    "b": str(b_value),
                    "cubic": str(cubic),
                    "cubic_discriminant": str(discriminant),
                    "square_discriminant": True,
                    "irreducible": True,
                }
                try:
                    branch_polynomial = (
                        discriminant.numerator() * discriminant.denominator()
                    )
                    collision = branch_polynomial.gcd(elliptic_bad)
                    hit["branch_collision_gcd"] = str(collision)
                    hit["branch_collision_detected"] = collision.degree() > 0
                except Exception as exc:
                    hit["branch_collision_error"] = repr(exc)
                    hit["branch_collision_detected"] = None
                try:
                    extension = function_field.extension(cubic, names=("z",))
                    hit["extension_genus"] = int(extension.genus())
                except Exception as exc:
                    hit["genus_error"] = repr(exc)
                hits.append(hit)
                report["hit_count"] += 1
                if len(hits) >= args.max_hits:
                    break
            trace_reports.append(report)
            if len(hits) >= args.max_hits:
                break

        output.update(
            {
                "status": "pass",
                "e8_root_count": len(roots),
                "e8_norm_six_trace_count": len(trace_vectors),
                "selected_trace_count": len(selected),
                "trace_reports": trace_reports,
                "hit_count": len(hits),
                "genus_one_smooth_branch_hit_count": sum(
                    hit.get("extension_genus") == 1
                    and hit.get("branch_collision_detected") is False
                    for hit in hits
                ),
                "hits": hits,
            }
        )
    except Exception as exc:
        output["status"] = "error"
        output["error"] = repr(exc)
        output["traceback"] = traceback.format_exc()

    output["elapsed_seconds"] = time.time() - started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": output["status"],
                "prime": args.prime,
                "trace_start": args.trace_start,
                "selected_traces": output.get("selected_trace_count"),
                "hits": output.get("hit_count"),
                "genus_one_smooth_hits": output.get(
                    "genus_one_smooth_branch_hit_count"
                ),
                "error": output.get("error"),
                "elapsed_seconds": output["elapsed_seconds"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
