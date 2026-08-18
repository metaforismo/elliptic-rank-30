#!/usr/bin/env python3
"""Fast replay wrapper for the three-variable additive-IV model.

The mathematical implementation lives in
rank17_iv_log_derivative_three_variable.py.  This wrapper replaces only the
lower-coefficient solver by the exact triangular formulas obtained from the
coefficients of E/P:

    m2 = 8*(-k+2*r3+4),

    m1 = -(k^2+14*k*r3+4*k+24*r2+13*r3^2-56*r3-80)/3,

    m0 = m0_base-16*r1.

The coefficients at degrees 10 and 9 are then a two-by-two affine system in
r0,r1.  Every resulting candidate is still passed through the full differential
identity, polynomial divisibility, A/B reconstruction, original six-parameter
reconstruction, and exact set comparison of the base module.

The leading recurrence coefficient is 28512=2^5*3^4*11, so F11 is a bad prime
for this coordinate chart even though the original surface locus is valid
there.  This wrapper accepts only primes outside {2,3,5,11}.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "rank17_iv_log_derivative_three_variable.py"
BAD_PRIMES = {2, 3, 5, 11}


def load_base():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_log_derivative_three_variable", BASE_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install_fast_solver(module):
    def solve_lower_variables(prime: int, r2: int, r3: int, k: int):
        if prime in BAD_PRIMES:
            return None, "BAD_PRIME_FOR_TRIANGULAR_RECURRENCE", 0

        inverse_three = pow(3, -1, prime)
        inverse_405 = pow(405 % prime, -1, prime)
        m2 = 8 * (-k + 2 * r3 + 4) % prime
        m1 = -(
            k * k
            + 14 * k * r3
            + 4 * k
            + 24 * r2
            + 13 * r3 * r3
            - 56 * r3
            - 80
        ) * inverse_three % prime
        numerator_without_r1 = (
            5 * k**3
            + 60 * k * k * r3
            + 30 * k * k
            + 1080 * k * r2
            - 120 * k * r3 * r3
            + 420 * k * r3
            + 456 * k
            + 2700 * r2 * r3
            - 4320 * r2
            - 715 * r3**3
            + 480 * r3 * r3
            - 6240 * r3
            - 9536
        )
        m0_base = -numerator_without_r1 * inverse_405 % prime

        def degree_10_9(r0_value: int, r1_value: int):
            m0_value = (m0_base - 16 * r1_value) % prime
            result = module.qbar_polynomial(
                prime=prime,
                r0=r0_value,
                r1=r1_value,
                r2=r2,
                r3=r3,
                k=k,
                m0=m0_value,
                m1=m1,
                m2=m2,
            )
            if result is None:
                raise AssertionError("E lost its exact P divisor")
            quotient = module.pad(result[0], 14, prime)
            return quotient[10] % prime, quotient[9] % prime

        constant_10, constant_9 = degree_10_9(0, 0)
        r0_10, r0_9 = degree_10_9(1, 0)
        r1_10, r1_9 = degree_10_9(0, 1)
        a = (r0_10 - constant_10) % prime
        c = (r0_9 - constant_9) % prime
        b = (r1_10 - constant_10) % prime
        d = (r1_9 - constant_9) % prime
        determinant = (a * d - b * c) % prime
        if determinant == 0:
            return None, "SINGULAR_THREE_BY_THREE_DETERMINANT", 0
        inverse_determinant = pow(determinant, -1, prime)
        rhs0 = -constant_10 % prime
        rhs1 = -constant_9 % prime
        r0 = (rhs0 * d - b * rhs1) * inverse_determinant % prime
        r1 = (a * rhs1 - rhs0 * c) * inverse_determinant % prime
        m0 = (m0_base - 16 * r1) % prime

        values = {
            "r0": r0,
            "r1": r1,
            "r2": r2,
            "r3": r3,
            "k": k,
            "m0": m0,
            "m1": m1,
            "m2": m2,
        }
        result = module.qbar_polynomial(prime=prime, **values)
        if result is None or any(result[0]):
            return (
                None,
                "NONZERO_REMAINING_DIFFERENTIAL_COEFFICIENT",
                determinant,
            )
        return values, "DIFFERENTIAL_EQUATION_ZERO", determinant

    module.solve_lower_variables = solve_lower_variables


def build(paths: list[Path]):
    module = load_base()
    install_fast_solver(module)
    payload = module.build(paths)
    payload["solver_variant"] = "exact triangular fast replay"
    payload["leading_recurrence_coefficient"] = 28512
    payload["leading_recurrence_factorization"] = "2^5*3^4*11"
    payload["bad_primes_for_this_coordinate_chart"] = sorted(BAD_PRIMES)
    payload["coordinate_chart_prime_boundary"] = (
        "F11 remains a valid prime for the original surface locus, but the "
        "degree-13 recurrence is inseparable there because 28512=0 mod 11."
    )
    payload.pop("record_sha256", None)
    payload["record_sha256"] = module.canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.census)
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
        "primes": [item["prime"] for item in payload["analyses"]],
        "surface_counts": {
            str(item["prime"]): item["source_surface_count"]
            for item in payload["analyses"]
        },
        "determinant_zero_source_counts": {
            str(item["prime"]): item["source_determinant_zero_count"]
            for item in payload["analyses"]
        },
        "bad_primes": payload["bad_primes_for_this_coordinate_chart"],
        "all_exact_matches": all(
            item["exact_surface_set_match"] for item in payload["analyses"]
        ),
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
