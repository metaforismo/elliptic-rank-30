#!/usr/bin/env python3
"""Three-variable logarithmic-derivative model of the additive-IV surface locus.

Put s=t-1 and, for a normalized surface,

    c4=s^2 A,  c6=s^2 B,
    H=s^2 A^3-B^2=kappa*t^3*P,

where A,B are monic of degrees 6,10 and P=t*R is monic of degree 5 with
R=t^4+r3*t^3+r2*t^2+r1*t+r0.  Define

    Q=P'+3*P/t.

From Q*H-P*H'=0 and gcd(s*A,B)=1 one gets a unique linear polynomial

    K=-12*t+k

such that

    B*K = A*(s*Q-2*P)-3*P*s*A',
    s*A^2*K = Q*B-2*P*B'.

Modulo P this implies A*K^2=Q^2+M*P with

    M=80*t^3+m2*t^2+m1*t+m0.

Writing N=Q^2+M*P and

    T=N*K*(s*Q-2*P)-3*P*s*(N'*K-2*N*K'),

the second equation becomes

    E=Q*T*K-2*P*T'*K+8*P*T*K'-s*N^2*K^2=0.

E is divisible by P.  The coefficients of E/P at degrees 13 and 12 determine
m2 and m1.  Degrees 11,10,9 are affine-linear in m0,r1,r0; away from their
explicit determinant divisor they determine those three values.  Thus the
remaining exact test is an exhaustion of only (r2,r3,k) in F_p^3.

This script derives the reduced coordinates from every certified finite-field
surface, exhausts the nonsingular three-variable chart, reconstructs A,B and
the original six parameters, and requires exact equality of the resulting
surface set with the source census.  Equality at finitely many primes is not a
characteristic-zero proof, and the singular determinant chart remains a
separate algebraic component to analyse.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def mod(value: int, prime: int) -> int:
    return value % prime


def trim(values: list[int], prime: int) -> list[int]:
    result = [value % prime for value in values]
    while result and result[-1] == 0:
        result.pop()
    return result


def pad(values: list[int], length: int, prime: int) -> list[int]:
    result = [value % prime for value in values]
    return result + [0] * max(0, length-len(result))


def add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    return trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ], prime)


def scale(values: list[int], scalar: int, prime: int) -> list[int]:
    return trim([scalar*value for value in values], prime)


def subtract(left: list[int], right: list[int], prime: int) -> list[int]:
    return add(left, scale(right, -1, prime), prime)


def multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    if not left or not right:
        return []
    result = [0] * (len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] = (result[i+j] + a*b) % prime
    return trim(result, prime)


def power(values: list[int], exponent: int, prime: int) -> list[int]:
    result = [1]
    base = list(values)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base, prime)
        base = multiply(base, base, prime)
        remaining //= 2
    return result


def derivative(values: list[int], prime: int) -> list[int]:
    return trim([
        index*values[index] for index in range(1, len(values))
    ], prime)


def divide_exact(
    dividend: list[int], divisor: list[int], prime: int
) -> list[int] | None:
    numerator = trim(dividend, prime)
    denominator = trim(divisor, prime)
    if not denominator:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(numerator) < len(denominator):
        return [] if not numerator else None
    quotient = [0] * (len(numerator)-len(denominator)+1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while numerator and len(numerator) >= len(denominator):
        shift = len(numerator)-len(denominator)
        coefficient = numerator[-1]*inverse_lead % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index+shift] = (
                numerator[index+shift]-coefficient*value
            ) % prime
        numerator = trim(numerator, prime)
    return trim(quotient, prime) if not numerator else None


def evaluate(values: list[int], argument: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(values):
        result = (result*argument+coefficient) % prime
    return result


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def coefficient(values: list[int], degree: int, prime: int) -> int:
    return values[degree] % prime if degree < len(values) else 0


def qbar_polynomial(
    *,
    prime: int,
    r0: int,
    r1: int,
    r2: int,
    r3: int,
    k: int,
    m0: int,
    m1: int,
    m2: int,
) -> tuple[list[int], dict[str, list[int]]] | None:
    s = [-1, 1]
    R = [r0, r1, r2, r3, 1]
    P = [0] + R
    Q = add(derivative(P, prime), scale(R, 3, prime), prime)
    K = [k, -12]
    Kprime = [-12 % prime]
    M = [m0, m1, m2, 80]
    N = add(power(Q, 2, prime), multiply(P, M, prime), prime)
    first = multiply(
        multiply(N, K, prime),
        subtract(multiply(s, Q, prime), scale(P, 2, prime), prime),
        prime,
    )
    correction = scale(
        multiply(
            multiply(P, s, prime),
            subtract(
                multiply(derivative(N, prime), K, prime),
                scale(multiply(N, Kprime, prime), 2, prime),
                prime,
            ),
            prime,
        ),
        3,
        prime,
    )
    T = subtract(first, correction, prime)
    E = subtract(
        add(
            subtract(
                multiply(multiply(Q, T, prime), K, prime),
                scale(
                    multiply(
                        multiply(P, derivative(T, prime), prime),
                        K,
                        prime,
                    ),
                    2,
                    prime,
                ),
                prime,
            ),
            scale(
                multiply(
                    multiply(P, T, prime), Kprime, prime
                ),
                8,
                prime,
            ),
            prime,
        ),
        multiply(
            multiply(s, power(N, 2, prime), prime),
            power(K, 2, prime),
            prime,
        ),
        prime,
    )
    quotient = divide_exact(E, P, prime)
    if quotient is None:
        return None
    return quotient, {
        "R": R,
        "P": P,
        "Q": Q,
        "K": K,
        "M": M,
        "N": N,
        "T": T,
        "E": E,
    }


def solve_affine_scalar(function, prime: int) -> int | None:
    f0 = function(0) % prime
    f1 = function(1) % prime
    slope = (f1-f0) % prime
    if function(2) % prime != (f0+2*slope) % prime:
        raise AssertionError("purported scalar equation is not affine-linear")
    if slope == 0:
        return None
    value = -f0*pow(slope, -1, prime) % prime
    if function(value) % prime:
        raise AssertionError("affine scalar solve did not kill the equation")
    return value


def solve_square_matrix(
    matrix: list[list[int]], rhs: list[int], prime: int
) -> tuple[list[int] | None, int]:
    size = len(matrix)
    work = [
        [value % prime for value in matrix[row]] + [rhs[row] % prime]
        for row in range(size)
    ]
    determinant = 1
    sign = 1
    for column in range(size):
        selected = next((
            row for row in range(column, size)
            if work[row][column]
        ), None)
        if selected is None:
            return None, 0
        if selected != column:
            work[column], work[selected] = work[selected], work[column]
            sign = -sign
        pivot = work[column][column]
        determinant = determinant*pivot % prime
        inverse = pow(pivot, -1, prime)
        work[column] = [value*inverse % prime for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (value-factor*pivot_value) % prime
                for value, pivot_value in zip(
                    work[row], work[column], strict=True
                )
            ]
    determinant = determinant*sign % prime
    return [work[row][-1] for row in range(size)], determinant


def solve_lower_variables(
    prime: int, r2: int, r3: int, k: int
) -> tuple[dict[str, int] | None, str, int]:
    def qbar(
        *, r0=0, r1=0, m0=0, m1=0, m2=0
    ) -> list[int]:
        result = qbar_polynomial(
            prime=prime,
            r0=r0,
            r1=r1,
            r2=r2,
            r3=r3,
            k=k,
            m0=m0,
            m1=m1,
            m2=m2,
        )
        if result is None:
            raise AssertionError("E lost its forced P factor")
        return pad(result[0], 14, prime)

    m2 = solve_affine_scalar(
        lambda value: coefficient(qbar(m2=value), 13, prime),
        prime,
    )
    if m2 is None:
        return None, "SINGULAR_M2_RECURRENCE", 0
    m1 = solve_affine_scalar(
        lambda value: coefficient(qbar(m2=m2, m1=value), 12, prime),
        prime,
    )
    if m1 is None:
        return None, "SINGULAR_M1_RECURRENCE", 0

    degrees = [11, 10, 9]
    base_values = [
        coefficient(qbar(m2=m2, m1=m1), degree, prime)
        for degree in degrees
    ]
    columns = []
    for name in ("m0", "r1", "r0"):
        kwargs = {"m2": m2, "m1": m1, name: 1}
        values = [
            coefficient(qbar(**kwargs), degree, prime)
            for degree in degrees
        ]
        columns.append([
            (value-base) % prime
            for value, base in zip(values, base_values, strict=True)
        ])
    matrix = [
        [columns[column][row] for column in range(3)]
        for row in range(3)
    ]

    for first in range(3):
        for second in range(first, 3):
            kwargs = {"m2": m2, "m1": m1}
            names = ("m0", "r1", "r0")
            kwargs[names[first]] = kwargs.get(names[first], 0)+1
            kwargs[names[second]] = kwargs.get(names[second], 0)+1
            values = [
                coefficient(qbar(**kwargs), degree, prime)
                for degree in degrees
            ]
            expected = [
                (
                    base_values[row]
                    + columns[first][row]
                    + columns[second][row]
                ) % prime
                for row in range(3)
            ]
            if values != expected:
                raise AssertionError("three-variable lower system is not affine")

    solution, determinant = solve_square_matrix(
        matrix,
        [(-value) % prime for value in base_values],
        prime,
    )
    if solution is None:
        return None, "SINGULAR_THREE_BY_THREE_DETERMINANT", 0
    m0, r1, r0 = solution
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
    quotient, _pieces = qbar_polynomial(prime=prime, **values)
    if quotient is None or any(quotient):
        return None, "NONZERO_REMAINING_DIFFERENTIAL_COEFFICIENT", determinant
    return values, "DIFFERENTIAL_EQUATION_ZERO", determinant


def build_original_c4(parameters: dict[str, int], prime: int) -> list[int]:
    p0 = parameters["p0"]
    p1 = parameters["p1"]
    p2 = parameters["p2"]
    p3 = parameters["p3"]
    r = parameters["r"]
    s = parameters["s"]
    values = [0] * 9
    values[0] = p0*p0 % prime
    values[1] = 2*p0*p1 % prime
    values[2] = (2*p0*p2+p1*p1) % prime
    values[3] = (2*p0*p3+2*p1*p2) % prime
    l0 = sum((values[index] for index in range(4)), 0)+s+1
    l1 = values[1]+2*values[2]+3*values[3]+7*s+8
    values[4] = (l1-5*l0+r) % prime
    values[5] = (4*l0-l1-2*r) % prime
    values[6] = r % prime
    values[7] = s % prime
    values[8] = 1
    return [value % prime for value in values]


def reconstruct_surface(
    prime: int, values: dict[str, int]
) -> dict[str, object] | None:
    result = qbar_polynomial(prime=prime, **values)
    if result is None:
        return None
    quotient, pieces = result
    if any(quotient):
        return None
    K2 = power(pieces["K"], 2, prime)
    K4 = power(pieces["K"], 4, prime)
    A = divide_exact(pieces["N"], K2, prime)
    B = divide_exact(pieces["T"], K4, prime)
    if A is None or B is None:
        return None
    if len(A) != 7 or A[-1] != 1 or len(B) != 11 or B[-1] != 1:
        return None

    s_polynomial = [-1, 1]
    X = subtract(
        multiply(
            A,
            subtract(
                multiply(s_polynomial, pieces["Q"], prime),
                scale(pieces["P"], 2, prime),
                prime,
            ),
            prime,
        ),
        scale(
            multiply(
                multiply(pieces["P"], s_polynomial, prime),
                derivative(A, prime),
                prime,
            ),
            3,
            prime,
        ),
        prime,
    )
    if X != multiply(B, pieces["K"], prime):
        return None
    Y = subtract(
        multiply(pieces["Q"], B, prime),
        scale(multiply(pieces["P"], derivative(B, prime), prime), 2, prime),
        prime,
    )
    expected_Y = multiply(
        multiply(s_polynomial, power(A, 2, prime), prime),
        pieces["K"],
        prime,
    )
    if Y != expected_Y:
        return None

    H = subtract(
        multiply(power(s_polynomial, 2, prime), power(A, 3, prime), prime),
        power(B, 2, prime),
        prime,
    )
    target = [0, 0, 0] + pieces["P"]
    if len(H) != len(target):
        return None
    scalar = H[-1] % prime
    if scalar == 0 or H != scale(target, scalar, prime):
        return None

    a0 = A[0] % prime
    b0 = B[0] % prime
    if a0 == 0:
        return None
    p0 = b0*pow(a0, -1, prime) % prime
    if p0 == 0 or p0*p0 % prime != a0 or p0*p0*p0 % prime != b0:
        return None
    inverse_2p0 = pow(2*p0 % prime, -1, prime)
    p1 = (A[1]-2*p0*p0)*inverse_2p0 % prime
    p2 = (
        A[2]-3*p0*p0-4*p0*p1-p1*p1
    )*inverse_2p0 % prime
    p3 = (
        A[3]-4*p0*p0-6*p0*p1-4*p0*p2-2*p1*p1-2*p1*p2
    )*inverse_2p0 % prime
    surface_r = (A[4]-2*A[5]+1) % prime
    surface_s = (A[5]-2) % prime
    parameters = {
        "p0": p0,
        "p1": p1,
        "p2": p2,
        "p3": p3,
        "r": surface_r,
        "s": surface_s,
    }
    c4 = multiply(power(s_polynomial, 2, prime), A, prime)
    c6 = multiply(power(s_polynomial, 2, prime), B, prime)
    if pad(c4, 9, prime) != build_original_c4(parameters, prime):
        return None
    if coefficient(c6, 0, prime) != p0**3 % prime:
        return None
    expected_jets = [
        3*p0*p0*p1,
        3*(p0*p0*p2+p0*p1*p1),
        3*p0*p0*p3+6*p0*p1*p2+p1**3,
    ]
    if any(
        coefficient(c6, index, prime) != expected_jets[index-1] % prime
        for index in range(1, 4)
    ):
        return None
    if evaluate(c6, 1, prime) or evaluate(derivative(c6, prime), 1, prime):
        return None

    record = {
        "reduced_coordinates": values,
        "linear_system_determinant": (
            -2*values["k"]**2
            -37*values["k"]*values["r3"]
            +52*values["k"]
            +108*values["r2"]
            -116*values["r3"]**2
            +148*values["r3"]
            -320
        ) % prime,
        "parameters": parameters,
        "A_coefficients_ascending": A,
        "B_coefficients_ascending": B,
        "residual_quartic_monic_coefficients_ascending": pieces["R"],
        "K_coefficients_ascending": pieces["K"],
        "M_coefficients_ascending": pieces["M"],
        "discriminant_scalar": scalar,
        "c4_coefficients_ascending": pad(c4, 9, prime),
        "c6_coefficients_ascending": pad(c6, 13, prime),
    }
    record["record_sha256"] = canonical_hash(record)
    return record


def derive_from_surface(record: dict[str, object], prime: int) -> dict[str, object]:
    c4 = [int(value) % prime for value in record["c4_coefficients_ascending"]]
    c6 = [int(value) % prime for value in record["c6_coefficients_ascending"]]
    s2 = [1, -2, 1]
    A = divide_exact(c4, s2, prime)
    B = divide_exact(c6, s2, prime)
    if A is None or B is None:
        raise AssertionError("surface invariants lost the IV square factor")
    H = subtract(
        multiply(s2, power(A, 3, prime), prime),
        power(B, 2, prime),
        prime,
    )
    if any(coefficient(H, degree, prime) for degree in range(4)):
        raise AssertionError("surface discriminant lost t^4")
    residual = trim(H[4:], prime)
    if len(residual) != 5 or residual[-1] == 0:
        raise AssertionError("surface residual is not quartic")
    inverse_lead = pow(residual[-1], -1, prime)
    R = scale(residual, inverse_lead, prime)
    P = [0] + R
    Q = add(derivative(P, prime), scale(R, 3, prime), prime)
    s_polynomial = [-1, 1]
    X = subtract(
        multiply(
            A,
            subtract(
                multiply(s_polynomial, Q, prime),
                scale(P, 2, prime),
                prime,
            ),
            prime,
        ),
        scale(
            multiply(
                multiply(P, s_polynomial, prime),
                derivative(A, prime),
                prime,
            ),
            3,
            prime,
        ),
        prime,
    )
    K = divide_exact(X, B, prime)
    if K is None or len(K) != 2 or K[1] != (-12) % prime:
        raise AssertionError(("surface did not produce K=-12t+k", K))
    Y = subtract(
        multiply(Q, B, prime),
        scale(multiply(P, derivative(B, prime), prime), 2, prime),
        prime,
    )
    if Y != multiply(
        multiply(s_polynomial, power(A, 2, prime), prime),
        K,
        prime,
    ):
        raise AssertionError("second logarithmic-derivative syzygy failed")
    M = divide_exact(
        subtract(multiply(A, power(K, 2, prime), prime), power(Q, 2, prime), prime),
        P,
        prime,
    )
    if M is None or len(M) != 4 or M[3] != 80 % prime:
        raise AssertionError(("surface did not produce cubic M with lead 80", M))
    values = {
        "r0": R[0],
        "r1": R[1],
        "r2": R[2],
        "r3": R[3],
        "k": K[0],
        "m0": M[0],
        "m1": M[1],
        "m2": M[2],
    }
    rebuilt = reconstruct_surface(prime, values)
    if rebuilt is None:
        raise AssertionError("derived reduced coordinates failed reconstruction")
    source_parameters = {
        key: int(value) % prime
        for key, value in record["parameters"].items()
        if key in {"p0", "p1", "p2", "p3", "r", "s"}
    }
    if rebuilt["parameters"] != source_parameters:
        raise AssertionError((
            "reduced reconstruction changed surface parameters",
            source_parameters,
            rebuilt["parameters"],
        ))
    return rebuilt


def parameter_tuple(record: dict[str, object]) -> tuple[int, ...]:
    parameters = record["parameters"]
    return tuple(
        int(parameters[name])
        for name in ("p0", "p1", "p2", "p3", "r", "s")
    )


def analyse_census(path: Path) -> dict[str, object]:
    source = json.loads(path.read_text(encoding="utf-8"))
    prime = int(source["prime"])
    if prime in (2, 3, 5):
        raise ValueError("the reduced recurrence uses denominators divisible by 2,3,5")
    source_records = source["records"]
    derived = [derive_from_surface(record, prime) for record in source_records]
    source_set = {parameter_tuple(record) for record in source_records}
    derived_set = {parameter_tuple(record) for record in derived}
    if source_set != derived_set or len(derived_set) != len(derived):
        raise AssertionError("surface-to-reduced map is not injective on the source census")

    determinant_zero_source_count = sum(
        record["linear_system_determinant"] == 0 for record in derived
    )
    reconstructed = []
    classifications: dict[str, int] = {}
    determinant_zero_key_count = 0
    for r2 in range(prime):
        for r3 in range(prime):
            for k in range(prime):
                values, classification, determinant = solve_lower_variables(
                    prime, r2, r3, k
                )
                classifications[classification] = (
                    classifications.get(classification, 0)+1
                )
                if determinant == 0:
                    determinant_zero_key_count += 1
                if values is None:
                    continue
                record = reconstruct_surface(prime, values)
                if record is not None:
                    reconstructed.append(record)

    reconstructed_set = {parameter_tuple(record) for record in reconstructed}
    duplicates = len(reconstructed)-len(reconstructed_set)
    if duplicates:
        raise AssertionError(("reduced chart reconstructs duplicate surfaces", duplicates))
    exact_match = reconstructed_set == source_set
    if not exact_match:
        raise AssertionError({
            "missing_from_reduced_chart": sorted(source_set-reconstructed_set),
            "extra_from_reduced_chart": sorted(reconstructed_set-source_set),
        })

    result = {
        "prime": prime,
        "source_record_sha256": source["record_sha256"],
        "source_surface_count": len(source_records),
        "derived_reduced_record_count": len(derived),
        "source_determinant_zero_count": determinant_zero_source_count,
        "three_variable_keys_visited": prime**3,
        "determinant_zero_key_count": determinant_zero_key_count,
        "classification_counts": classifications,
        "reconstructed_surface_count": len(reconstructed),
        "exact_surface_set_match": True,
        "derived_records": derived,
        "reconstructed_record_sha256": [
            record["record_sha256"] for record in reconstructed
        ],
    }
    result["record_sha256"] = canonical_hash(result)
    return result


def build(paths: list[Path]) -> dict[str, object]:
    analyses = [analyse_census(path) for path in paths]
    analyses.sort(key=lambda item: item["prime"])
    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_log_derivative_three_variable",
        "truth_status": (
            "EXACT finite-field equivalence between the certified full surface "
            "censuses and the nonsingular three-variable logarithmic-derivative "
            "chart at the supplied primes; the singular determinant divisor and "
            "characteristic-zero geometry remain unresolved; no section or rank-30 conclusion"
        ),
        "reduced_variables": ["r2", "r3", "k"],
        "reconstructed_variables": ["m2", "m1", "m0", "r1", "r0"],
        "linear_polynomial": "K(t)=-12*t+k",
        "cubic_polynomial": "M(t)=80*t^3+m2*t^2+m1*t+m0",
        "differential_equation": (
            "With N=Q^2+M*P and T=N*K*((t-1)*Q-2*P)-3*P*(t-1)*(N'*K-2*N*K'), "
            "the exact equation is Q*T*K-2*P*T'*K+8*P*T*K'-(t-1)*N^2*K^2=0."
        ),
        "analyses": analyses,
        "limitations": [
            "The three-by-three determinant divisor is not eliminated as a characteristic-zero component.",
            "Finite-field equality at finitely many primes does not prove equality of schemes over Q.",
            "The residual squarefree and split open conditions are not imposed because the target here is the full surface curve.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.census)
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True)+"\n",
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
        "all_exact_matches": all(
            item["exact_surface_set_match"] for item in payload["analyses"]
        ),
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
