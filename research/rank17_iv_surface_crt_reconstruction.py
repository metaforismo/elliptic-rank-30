#!/usr/bin/env python3
"""Bounded exact CRT reconstruction of normalized split IV surfaces.

Pairs every geometric representative over F11 and F13.  For each coordinate it
tries the balanced integer modulo 143 and every reduced rational a/b with
|a|, b <= the declared bound and b invertible modulo 143.  Candidates are
promoted only after a complete exact Q verification of the five surface
equations, exact I4/IV/I12 orders, split tangent squares, and squarefree residual
discriminant.

Failure to reconstruct inside the stated bound is not a global obstruction.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp

PARAMETER_NAMES = ("p0", "p1", "p2", "p3", "r", "s")


def convolution(left, right, length=None):
    if length is None:
        length = len(left) + len(right) - 1
    zero = left[0] * 0
    result = [zero for _ in range(length)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i + j >= length:
                break
            result[i + j] += left_value * right_value
    return result


def trim(values):
    result = list(values)
    while result and result[-1] == 0:
        result.pop()
    return result


def polynomial_subtract(left, right):
    size = max(len(left), len(right))
    return trim([
        (left[index] if index < len(left) else 0)
        - (right[index] if index < len(right) else 0)
        for index in range(size)
    ])


def polynomial_power(values, exponent):
    result = [Fraction(1)]
    base = list(values)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = convolution(result, base)
        remaining //= 2
        if remaining:
            base = convolution(base, base)
    return trim(result)


def evaluate(values, argument):
    result = Fraction(0)
    for coefficient in reversed(values):
        result = result * argument + coefficient
    return result


def derivative_at_one(values):
    return sum(
        (Fraction(index) * values[index] for index in range(1, len(values))),
        Fraction(0),
    )


def quadratic_coefficient_at_one(values):
    return sum(
        (
            Fraction(index * (index - 1), 2) * values[index]
            for index in range(2, len(values))
        ),
        Fraction(0),
    )


def build_invariants(parameters: dict[str, Fraction]):
    p0 = parameters["p0"]
    p1 = parameters["p1"]
    p2 = parameters["p2"]
    p3 = parameters["p3"]
    r = parameters["r"]
    s = parameters["s"]

    a0 = p0**2
    a1 = 2*p0*p1
    a2 = 2*p0*p2 + p1**2
    a3 = 2*p0*p3 + 2*p1*p2
    l0 = a0 + a1 + a2 + a3 + s + 1
    l1 = a1 + 2*a2 + 3*a3 + 7*s + 8
    a4 = l1 - 5*l0 + r
    a5 = 4*l0 - l1 - 2*r
    c4 = [a0, a1, a2, a3, a4, a5, r, s, Fraction(1)]

    reversed_c4 = list(reversed(c4))
    cube = convolution(
        convolution(reversed_c4, reversed_c4, 12),
        reversed_c4,
        12,
    )
    root = [Fraction(0) for _ in range(12)]
    root[0] = Fraction(1)
    for order in range(1, 12):
        correction = sum(
            (root[index] * root[order-index] for index in range(1, order)),
            Fraction(0),
        )
        root[order] = (cube[order] - correction) / 2

    c6 = [Fraction(0) for _ in range(13)]
    c6[0] = p0**3
    c6[12] = Fraction(1)
    for order in range(1, 12):
        c6[12-order] = root[order]

    equations = [
        c6[1] - 3*p0**2*p1,
        c6[2] - 3*(p0**2*p2 + p0*p1**2),
        c6[3] - (3*p0**2*p3 + 6*p0*p1*p2 + p1**3),
        evaluate(c6, 1),
        derivative_at_one(c6),
    ]
    delta = polynomial_subtract(
        polynomial_power(c4, 3),
        polynomial_power(c6, 2),
    )
    delta += [Fraction(0)] * (25 - len(delta))
    return c4, c6, equations, delta


def rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = value.numerator
    denominator = value.denominator
    numerator_root = math.isqrt(numerator)
    denominator_root = math.isqrt(denominator)
    if numerator_root**2 != numerator or denominator_root**2 != denominator:
        return None
    return Fraction(numerator_root, denominator_root)


def sympy_poly(values, symbol):
    expression = sum(
        sp.Rational(value.numerator, value.denominator) * symbol**index
        for index, value in enumerate(values)
    )
    return sp.Poly(expression, symbol, domain=sp.QQ)


def exact_divide_by_fixed_factors(delta):
    t = sp.symbols("t")
    delta_poly = sympy_poly(delta, t)
    fixed = sp.Poly(t**4 * (t-1)**4, t, domain=sp.QQ)
    quotient, remainder = sp.div(delta_poly, fixed, domain=sp.QQ)
    if not remainder.is_zero:
        raise AssertionError("discriminant is not divisible by t^4(t-1)^4")
    return t, quotient


def verify_over_q(parameters: dict[str, Fraction]) -> dict[str, object] | None:
    if parameters["p0"] == 0:
        return None
    c4, c6, equations, delta = build_invariants(parameters)
    if any(equations):
        return None
    if any(delta[index] for index in range(4)) or delta[4] == 0:
        return None
    if any(delta[index] for index in range(13, 25)) or delta[12] == 0:
        return None
    if evaluate(c4, 1) != 0 or derivative_at_one(c4) != 0:
        return None
    if evaluate(c6, 1) != 0 or derivative_at_one(c6) != 0:
        return None
    b_at_one = quadratic_coefficient_at_one(c6)
    if b_at_one == 0:
        return None

    i4_root = rational_square_root(-3 * parameters["p0"])
    iv_root = rational_square_root(-2 * b_at_one)
    if i4_root is None or iv_root is None:
        return None

    t, residual = exact_divide_by_fixed_factors(delta)
    if residual.degree() != 4:
        return None
    if sp.gcd(residual, residual.diff()).degree() != 0:
        return None
    c4_poly = sympy_poly(c4, t)
    if sp.gcd(c4_poly, residual).degree() != 0:
        return None

    record = {
        "parameters": {
            name: str(parameters[name]) for name in PARAMETER_NAMES
        },
        "c4_coefficients_ascending": [str(value) for value in c4],
        "c6_coefficients_ascending": [str(value) for value in c6],
        "five_surface_equations": [str(value) for value in equations],
        "discriminant": {
            "factorization": str(sympy_poly(delta, t).factor_list()),
            "residual_quartic": str(residual.as_expr()),
            "residual_squarefree": True,
            "residual_coprime_to_c4": True,
        },
        "split_tangents": {
            "i4_square_target": str(-3 * parameters["p0"]),
            "i4_root": str(i4_root),
            "iv_square_target": str(-2 * b_at_one),
            "iv_root": str(iv_root),
        },
    }
    raw = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    record["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return record


def crt_pair(left: int, right: int, left_modulus: int, right_modulus: int) -> int:
    modulus = left_modulus * right_modulus
    step = ((right - left) * pow(left_modulus, -1, right_modulus)) % right_modulus
    return (left + left_modulus * step) % modulus


def balanced(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if residue > modulus // 2 else residue


def bounded_rational_options(residue: int, modulus: int, bound: int):
    options = set()
    for denominator in range(1, bound + 1):
        if math.gcd(denominator, modulus) != 1:
            continue
        numerator = balanced(residue * denominator, modulus)
        if abs(numerator) > bound:
            continue
        if math.gcd(numerator, denominator) != 1:
            continue
        options.add(Fraction(numerator, denominator))
    return sorted(options)


def representative_residues(candidate):
    params = candidate["parameters"]
    return {name: int(params[name]) for name in PARAMETER_NAMES}


def build(
    f11_path: Path,
    f13_path: Path,
    rational_bound: int,
    combination_cap: int,
):
    f11 = json.loads(f11_path.read_text(encoding="utf-8"))
    f13 = json.loads(f13_path.read_text(encoding="utf-8"))
    if int(f11["prime"]) != 11 or int(f13["prime"]) != 13:
        raise ValueError("expected the F11 and F13 geometric quotient certificates")
    modulus = 11 * 13
    exact_records: dict[str, dict[str, object]] = {}
    pair_records = []

    for index11, candidate11 in enumerate(f11["representatives"]):
        residues11 = representative_residues(candidate11)
        for index13, candidate13 in enumerate(f13["representatives"]):
            residues13 = representative_residues(candidate13)
            crt = {
                name: crt_pair(residues11[name], residues13[name], 11, 13)
                for name in PARAMETER_NAMES
            }
            balanced_parameters = {
                name: Fraction(balanced(crt[name], modulus))
                for name in PARAMETER_NAMES
            }
            tested = 0
            promoted = []
            direct = verify_over_q(balanced_parameters)
            tested += 1
            if direct is not None:
                exact_records[direct["record_sha256"]] = direct
                promoted.append(direct["record_sha256"])

            options = {
                name: bounded_rational_options(
                    crt[name], modulus, rational_bound
                )
                for name in PARAMETER_NAMES
            }
            option_count = math.prod(len(options[name]) for name in PARAMETER_NAMES)
            truncated = option_count > combination_cap
            if not truncated:
                for values in itertools.product(
                    *(options[name] for name in PARAMETER_NAMES)
                ):
                    parameters = dict(zip(PARAMETER_NAMES, values, strict=True))
                    if parameters == balanced_parameters:
                        continue
                    tested += 1
                    record = verify_over_q(parameters)
                    if record is not None:
                        exact_records[record["record_sha256"]] = record
                        promoted.append(record["record_sha256"])
            pair_records.append({
                "f11_representative_index": index11,
                "f13_representative_index": index13,
                "crt_residues_mod_143": crt,
                "balanced_integer_parameters": {
                    name: str(value) for name, value in balanced_parameters.items()
                },
                "bounded_option_counts": {
                    name: len(options[name]) for name in PARAMETER_NAMES
                },
                "cartesian_option_count": option_count,
                "cartesian_search_truncated": truncated,
                "exact_candidates_tested": tested,
                "promoted_record_sha256": sorted(set(promoted)),
            })

    exact = [exact_records[key] for key in sorted(exact_records)]
    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_surface_crt_f11_f13_bounded",
        "truth_status": (
            "EXACT bounded CRT reconstruction and Q substitution; absence outside "
            "the declared bounds or truncated pair products is unresolved; no section "
            "or rank-30 conclusion"
        ),
        "source_quotient_sha256": {
            "11": f11["record_sha256"],
            "13": f13["record_sha256"],
        },
        "modulus": modulus,
        "rational_numerator_denominator_bound": rational_bound,
        "per_pair_combination_cap": combination_cap,
        "f11_geometric_surface_count": f11["geometric_surface_count"],
        "f13_geometric_surface_count": f13["geometric_surface_count"],
        "surface_pair_count": len(pair_records),
        "exact_q_surface_count": len(exact),
        "exact_q_surfaces": exact,
        "pairs": pair_records,
        "limitations": [
            "Coordinatewise CRT does not supply a canonical correspondence between different finite-field points.",
            "Only balanced integers and rationals within the declared bound are tested.",
            "A rational IV surface without the height-79/12 section does not realize the rank-17 seed lattice.",
            "The global unconditional rank lower bound remains 29.",
        ],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f11", type=Path, required=True)
    parser.add_argument("--f13", type=Path, required=True)
    parser.add_argument("--bound", type=int, default=8)
    parser.add_argument("--combination-cap", type=int, default=200000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build(
        arguments.f11,
        arguments.f13,
        arguments.bound,
        arguments.combination_cap,
    )
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if committed != payload:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
    print(json.dumps({
        "surface_pair_count": payload["surface_pair_count"],
        "exact_q_surface_count": payload["exact_q_surface_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
