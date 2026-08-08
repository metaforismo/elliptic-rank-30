#!/usr/bin/env python3
"""Exact elementary certificate for the collision ninth-power descent.

A nontrivial rational collision in the degree-(2,4) parameter family yields a
primitive integral solution of

    (a^3+b^3)(a^6+b^6) = 2^7 3^5 c^9.

This verifier checks the algebraic gcd identities and the complete residue
calculations that force the four simultaneous ninth-power equations recorded in
the accompanying proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from research.degree24_polynomial_section_obstruction import _variables, canonical_payload


def verify_gcd_identities() -> None:
    a, b = _variables(2)
    A = a + b
    B = a**2 - a * b + b**2
    C = a**2 + b**2
    D = a**4 - a**2 * b**2 + b**4

    identities = [
        B - A * (a - 2 * b) - 3 * b**2,
        C - A * (a - b) - 2 * b**2,
        D - A * a**2 * (a - b) - b**4,
        C - B - a * b,
        D - (a**2 + a * b - b**2) * B - 2 * b**3 * (b - a),
        D - (a**2 - 2 * b**2) * C - 3 * b**4,
    ]
    if any(not identity.is_zero() for identity in identities):
        raise AssertionError("one of the exact gcd identities failed")

    if A * B != a**3 + b**3:
        raise AssertionError("cubic cyclotomic factorization failed")
    if C * D != a**6 + b**6:
        raise AssertionError("sextic cyclotomic factorization failed")


def verify_two_adic_residues() -> None:
    for a in range(1, 8, 2):
        for b in range(1, 8, 2):
            B = a * a - a * b + b * b
            C = a * a + b * b
            D = a**4 - a * a * b * b + b**4
            if B % 2 != 1 or D % 2 != 1:
                raise AssertionError("B or D was not odd")
            if C % 8 != 2:
                raise AssertionError("v_2(C) is not exactly one")


def verify_three_adic_residues() -> None:
    checked = 0
    for a in range(9):
        for b in range(9):
            if a % 3 == 0 or b % 3 == 0 or (a + b) % 3 != 0:
                continue
            checked += 1
            B = a * a - a * b + b * b
            C = a * a + b * b
            D = a**4 - a * a * b * b + b**4
            if B % 3 != 0 or B % 9 == 0:
                raise AssertionError("v_3(B) is not exactly one")
            if C % 3 == 0 or D % 3 == 0:
                raise AssertionError("C or D was unexpectedly divisible by 3")
    if checked != 18:
        raise AssertionError(f"unexpected residue-class count: {checked}")


def build_certificate() -> dict:
    verify_gcd_identities()
    verify_two_adic_residues()
    verify_three_adic_residues()
    return {
        "certificate_id": "degree24-collision-ninth-power-descent-v1",
        "claim_status": "proved",
        "collision_equation": "(a^3+b^3)(a^6+b^6)=2^7*3^5*c^9",
        "coprimality": {
            "A_B": "gcd(a+b,a^2-ab+b^2) divides 3",
            "A_C": "gcd(a+b,a^2+b^2) divides 2",
            "A_D": "gcd(a+b,a^4-a^2b^2+b^4)=1",
            "B_C": "gcd(a^2-ab+b^2,a^2+b^2)=1",
            "B_D": "gcd(a^2-ab+b^2,a^4-a^2b^2+b^4)=1",
            "C_D": "gcd(a^2+b^2,a^4-a^2b^2+b^4)=1",
            "scope": "for primitive a,b with both odd and a+b divisible by 3",
        },
        "cyclotomic_factors": {
            "A": "a+b",
            "B": "a^2-ab+b^2",
            "C": "a^2+b^2",
            "D": "a^4-a^2b^2+b^4",
        },
        "forced_equations": {
            "A": "a+b=2^6*3^4*w^9",
            "B": "a^2-ab+b^2=3*x^9",
            "C": "a^2+b^2=2*y^9",
            "D": "a^4-a^2b^2+b^4=z^9",
        },
        "interpretation": {
            "eisenstein": "B is the norm of a+b*omega in Z[omega]",
            "gaussian": "C is the norm of a+i*b in Z[i]",
            "research_reduction": "a collision requires simultaneous ninth-power norm conditions in both UFDs",
        },
        "primitive_local_conditions": {
            "mod_2": "a and b are odd; v_2(C)=1; B and D are odd",
            "mod_3": "3 does not divide ab, a+b is divisible by 3, v_3(B)=1, and 3 does not divide CD",
            "valuation_A": "v_2(A)=6 mod 9 and v_3(A)=4 mod 9",
        },
        "schema_version": 1,
        "scope_limitation": "The descent is a necessary condition for a collision; nonexistence of the simultaneous ninth-power solutions is not yet proved here.",
        "verification": {
            "gcd_identities": "exact symbolic polynomial arithmetic",
            "mod_2": "complete enumeration of odd residue classes modulo 8",
            "mod_3": "complete enumeration of admissible residue classes modulo 9",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()

    payload = canonical_payload(build_certificate())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check and args.check.read_text(encoding="utf-8") != payload:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print("VERIFIED degree-(2,4) collision ninth-power descent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
