#!/usr/bin/env python3
"""Certify the residual discriminant geometry of finite-field IV surfaces.

The surface-locus enumerator enforces exact discriminant orders 4 at t=0,
4 at t=1, and 12 at infinity.  To obtain precisely the intended

    I12 + I4 + IV + 4 I1

configuration, the residual quartic R in

    Delta = t^4 (t-1)^4 R(t)

must have degree four, be squarefree, and be coprime to c4.  This script checks
those conditions exactly over F_p for every geometric representative and
produces a filtered certificate.  It makes no lift, section, or rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def trim(values: list[int], prime: int) -> list[int]:
    result = [value % prime for value in values]
    while result and result[-1] == 0:
        result.pop()
    return result


def degree(values: list[int], prime: int) -> int:
    return len(trim(values, prime)) - 1


def derivative(values: list[int], prime: int) -> list[int]:
    return trim(
        [index * values[index] for index in range(1, len(values))],
        prime,
    )


def divmod_poly(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    numerator = trim(dividend, prime)
    denominator = trim(divisor, prime)
    if not denominator:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(numerator) < len(denominator):
        return [], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while numerator and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] = (
                numerator[index + shift] - coefficient * value
            ) % prime
        numerator = trim(numerator, prime)
    return trim(quotient, prime), numerator


def gcd_poly(left: list[int], right: list[int], prime: int) -> list[int]:
    a = trim(left, prime)
    b = trim(right, prime)
    while b:
        _quotient, remainder = divmod_poly(a, b, prime)
        a, b = b, remainder
    if not a:
        return []
    inverse = pow(a[-1], -1, prime)
    return trim([inverse * value for value in a], prime)


def multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    if not left or not right:
        return []
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return trim(result, prime)


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def analyze_candidate(candidate: dict[str, object], prime: int) -> dict[str, object]:
    delta = [
        int(value) % prime
        for value in candidate["discriminant_coefficients_ascending"]
    ]
    c4 = [
        int(value) % prime
        for value in candidate["c4_coefficients_ascending"]
    ]
    t4 = [0, 0, 0, 0, 1]
    t_minus_one_4 = [1, -4, 6, -4, 1]
    fixed = multiply(t4, t_minus_one_4, prime)
    residual, remainder = divmod_poly(delta, fixed, prime)
    residual = trim(residual, prime)
    residual_derivative_gcd = gcd_poly(
        residual, derivative(residual, prime), prime
    )
    residual_c4_gcd = gcd_poly(residual, c4, prime)
    checks = {
        "fixed_factor_divides": not remainder,
        "residual_degree": degree(residual, prime),
        "residual_degree_four": degree(residual, prime) == 4,
        "residual_squarefree": residual_derivative_gcd == [1],
        "residual_coprime_to_c4": residual_c4_gcd == [1],
    }
    accepted = all((
        checks["fixed_factor_divides"],
        checks["residual_degree_four"],
        checks["residual_squarefree"],
        checks["residual_coprime_to_c4"],
    ))
    record = {
        "representative_record_sha256": candidate["record_sha256"],
        "parameters": candidate["parameters"],
        "residual_quartic_coefficients_ascending": residual,
        "gcd_residual_derivative": residual_derivative_gcd,
        "gcd_residual_c4": residual_c4_gcd,
        "checks": checks,
        "accepted_exact_fibre_configuration": accepted,
    }
    record["record_sha256"] = canonical_hash(record)
    return record


def build(quotient_path: Path) -> dict[str, object]:
    quotient = json.loads(quotient_path.read_text(encoding="utf-8"))
    prime = int(quotient["prime"])
    analyses = [
        analyze_candidate(candidate, prime)
        for candidate in quotient["representatives"]
    ]
    accepted_hashes = {
        record["representative_record_sha256"]
        for record in analyses
        if record["accepted_exact_fibre_configuration"]
    }
    accepted_representatives = [
        candidate
        for candidate in quotient["representatives"]
        if candidate["record_sha256"] in accepted_hashes
    ]
    rejected = [
        record for record in analyses
        if not record["accepted_exact_fibre_configuration"]
    ]
    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": f"rank17_iv_surface_residual_geometry_f{prime}",
        "truth_status": (
            f"EXACT residual-discriminant filtering over F_{prime}; "
            "no characteristic-zero, section, or rank-30 conclusion"
        ),
        "prime": prime,
        "source_geometric_quotient_sha256": quotient["record_sha256"],
        "input_geometric_surface_count": quotient["geometric_surface_count"],
        "accepted_exact_configuration_count": len(accepted_representatives),
        "rejected_degenerate_count": len(rejected),
        "exact_configuration": "I12 + I4 + IV + 4 I1",
        "analyses": analyses,
        "representatives": accepted_representatives,
        "limitations": [
            "The certificate is confined to the normalized finite-field chart.",
            "A finite-field surface need not lift to characteristic zero.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quotient", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build(arguments.quotient)
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
        "prime": payload["prime"],
        "input_geometric_surface_count": payload["input_geometric_surface_count"],
        "accepted_exact_configuration_count": payload["accepted_exact_configuration_count"],
        "rejected_degenerate_count": payload["rejected_degenerate_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
