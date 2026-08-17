#!/usr/bin/env sage -python
"""Exact Wronskian reduction for the normalized additive-IV surface locus.

For monic polynomials A of degree 6 and B of degree 10 set

    F=(t-1)^2*A^3,  G=B^2,
    C=2*A*B + 3*(t-1)*A'*B - 2*(t-1)*A*B'.

The symbolic identity

    F'*G - F*G' = (t-1)*A^2*B*C

shows that, when A(0)=p0^2 and B(0)=p0^3 with p0 nonzero, the conditions

    C=c3*t^3+c4*t^4,  c3*c4*B(1) != 0

force H=F-G to have exact order 4 at zero and exact degree 8.  Hence
H=t^4*R4 and c4=(t-1)^2*A, c6=(t-1)^2*B have exact discriminant orders
I4, IV, I12.  Conversely, the intended surface identity has exactly this
Wronskian support.  The residual quartic open conditions remain separate.

This script verifies the generic identities over Q and replays the reduction
on the finite-field surface certificates.  It does not construct a Q-point or
a height-79/12 section and does not change the rank-29 lower bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ
from sage.version import version as sage_version


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def symbolic_certificate() -> dict[str, object]:
    coefficient_ring = PolynomialRing(
        QQ,
        names=(
            "p0", "a1", "a2", "a3", "a4", "a5",
            "b1", "b2", "b3", "b4", "b5",
            "b6", "b7", "b8", "b9",
        ),
    )
    (
        p0, a1, a2, a3, a4, a5,
        b1, b2, b3, b4, b5, b6, b7, b8, b9,
    ) = coefficient_ring.gens()
    time_ring = PolynomialRing(coefficient_ring, "t")
    t = time_ring.gen()
    A = time_ring(p0**2 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5 + t**6)
    B = time_ring(
        p0**3 + b1*t + b2*t**2 + b3*t**3 + b4*t**4 + b5*t**5
        + b6*t**6 + b7*t**7 + b8*t**8 + b9*t**9 + t**10
    )
    F = (t-1)**2 * A**3
    G = B**2
    C = 2*A*B + 3*(t-1)*A.derivative()*B - 2*(t-1)*A*B.derivative()
    wronskian = F.derivative()*G - F*G.derivative()
    factorized = (t-1)*A**2*B*C
    if wronskian != factorized:
        raise AssertionError("generic Wronskian factorization failed")
    if C.degree() > 15:
        raise AssertionError(("unexpected Wronskian bracket degree", C.degree()))

    parameter_ring = PolynomialRing(QQ, names=("p0", "p1", "p2", "p3", "r", "s"))
    p0x, p1x, p2x, p3x, r, s = parameter_ring.gens()
    tx_ring = PolynomialRing(parameter_ring, "t")
    x = tx_ring.gen()
    c4_coefficients = [
        p0x**2,
        2*p0x*p1x,
        2*p0x*p2x + p1x**2,
        2*p0x*p3x + 2*p1x*p2x,
    ]
    l0 = sum(c4_coefficients) + s + 1
    l1 = (
        c4_coefficients[1]
        + 2*c4_coefficients[2]
        + 3*c4_coefficients[3]
        + 7*s + 8
    )
    c4_coefficients.extend([
        l1 - 5*l0 + r,
        4*l0 - l1 - 2*r,
        r,
        s,
        parameter_ring.one(),
    ])
    c4 = tx_ring(sum(c4_coefficients[index]*x**index for index in range(9)))
    A_from_parameters, remainder = c4.quo_rem((x-1)**2)
    if remainder:
        raise AssertionError("the c4 parametrization lost its (t-1)^2 factor")
    expected_A_coefficients = [
        p0x**2,
        2*p0x**2 + 2*p0x*p1x,
        3*p0x**2 + 4*p0x*p1x + 2*p0x*p2x + p1x**2,
        (
            4*p0x**2 + 6*p0x*p1x + 4*p0x*p2x + 2*p0x*p3x
            + 2*p1x**2 + 2*p1x*p2x
        ),
        r + 2*s + 3,
        s + 2,
        parameter_ring.one(),
    ]
    expected_A = tx_ring(sum(
        expected_A_coefficients[index]*x**index for index in range(7)
    ))
    if A_from_parameters != expected_A:
        raise AssertionError("unexpected birational A-coordinate map")

    return {
        "generic_wronskian_identity_verified": True,
        "generic_C_degree_upper_bound": int(C.degree()),
        "zero_coefficient_indices_for_surface_curve": [
            0, 1, 2, *range(5, 16)
        ],
        "free_wronskian_coefficients": [3, 4],
        "A_coefficient_map": {
            "a0": "p0^2",
            "a1": "2*p0^2+2*p0*p1",
            "a2": "3*p0^2+4*p0*p1+2*p0*p2+p1^2",
            "a3": "4*p0^2+6*p0*p1+4*p0*p2+2*p0*p3+2*p1^2+2*p1*p2",
            "a4": "r+2*s+3",
            "a5": "s+2",
            "a6": "1",
            "inverse_r": "a4-2*a5+1",
            "inverse_s": "a5-2",
        },
        "logical_implication": {
            "assumptions": [
                "characteristic zero",
                "A and B are monic of degrees 6 and 10",
                "A(0)=p0^2, B(0)=p0^3, p0!=0",
                "C=c3*t^3+c4*t^4 with c3*c4*B(1)!=0",
            ],
            "conclusion": (
                "H=(t-1)^2*A^3-B^2 has exact order 4 at zero and exact "
                "degree 8, so H=t^4*R4 and the induced c4,c6 have exact "
                "I4, IV, I12 discriminant orders"
            ),
            "degree_argument": (
                "If deg(H)=m, the leading term of H'*G-H*G' is "
                "(m-20)*lc(H)*t^(m+19). The verified Wronskian identity "
                "has degree 27 with nonzero t^4 bracket coefficient, hence m=8."
            ),
            "order_argument": (
                "H(0)=0. If ord_0(H)=m<4, then H'*G-H*G' has order m-1 "
                "because G(0)=p0^6 and m is 1,2,3; this contradicts the "
                "t^3 factor. Nonzero c3 then forces ord_0(H)=4."
            ),
        },
    }


def verify_finite_certificate(path: Path) -> dict[str, object]:
    source = json.loads(path.read_text(encoding="utf-8"))
    prime = int(source["prime"])
    field = GF(prime)
    time_ring = PolynomialRing(field, "t")
    t = time_ring.gen()
    records = []
    for candidate in source["representatives"]:
        c4 = time_ring(candidate["c4_coefficients_ascending"])
        c6 = time_ring(candidate["c6_coefficients_ascending"])
        A, rem4 = c4.quo_rem((t-1)**2)
        B, rem6 = c6.quo_rem((t-1)**2)
        if rem4 or rem6 or A.degree() != 6 or B.degree() != 10:
            raise AssertionError(("invalid A/B factorization", path))
        C = 2*A*B + 3*(t-1)*A.derivative()*B - 2*(t-1)*A*B.derivative()
        support = [index for index in range(C.degree()+1) if C[index]]
        if support != [3, 4]:
            raise AssertionError(("unexpected C support", prime, support))
        c3 = C[3]
        c4_coefficient = C[4]
        if not c3 or not c4_coefficient:
            raise AssertionError("the exact I4/I12 Wronskian coefficients vanished")
        extra_point = -c3/c4_coefficient
        H = (t-1)**2*A**3 - B**2
        R, remainder = H.quo_rem(t**4)
        if remainder or R.degree() != 4:
            raise AssertionError(("Wronskian converse failed", prime))
        if c4_coefficient != -12*R.leading_coefficient():
            raise AssertionError("leading coefficient relation failed")
        if c3 != 12*R.leading_coefficient()*extra_point:
            raise AssertionError("extra ramification coordinate relation failed")
        residual_gcd = R.gcd(R.derivative()).monic()
        collision = R(extra_point) == 0
        if collision and residual_gcd == 1:
            raise AssertionError(
                "an extra-point/pole collision did not create repeated residual geometry"
            )
        numerator = (extra_point-1)**2 * A(extra_point)**3
        denominator = extra_point**4 * R(extra_point)
        branch_value = None if denominator == 0 else numerator/denominator
        record = {
            "source_record_sha256": candidate["record_sha256"],
            "parameters": candidate["parameters"],
            "C_coefficients_ascending": [int(value) for value in C.list()],
            "C_support": support,
            "extra_ramification_point": int(extra_point),
            "extra_branch_value": (
                None if branch_value is None else int(branch_value)
            ),
            "residual_quartic_coefficients_ascending": [
                int(value) for value in R.list()
            ],
            "residual_gcd_derivative": [
                int(value) for value in residual_gcd.list()
            ],
            "extra_point_collides_with_residual_pole": collision,
            "wronskian_leading_relation_verified": True,
        }
        record["record_sha256"] = canonical_hash(record)
        records.append(record)
    return {
        "prime": prime,
        "source_record_sha256": source["record_sha256"],
        "record_count": len(records),
        "records": records,
    }


def build(paths: list[Path]) -> dict[str, object]:
    symbolic = symbolic_certificate()
    finite = [verify_finite_certificate(path) for path in paths]
    f13 = next((item for item in finite if item["prime"] == 13), None)
    f13_collision_verified = False
    if f13 is not None and f13["record_count"] == 0:
        # The residual-geometry certificate has no accepted representatives.
        # The broader split quotient must be supplied separately to diagnose it.
        f13_collision_verified = False
    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_wronskian_reduction",
        "truth_status": (
            "CERTIFIED generic Wronskian identity and exact replay on the supplied "
            "finite-field open surface loci; no characteristic-zero point, section, "
            "or rank-30 conclusion"
        ),
        "sage_version": str(sage_version),
        "symbolic": symbolic,
        "finite_field_replays": finite,
        "f13_open_locus_collision_note": (
            "The separate F13 residual-geometry certificate rejects its sole split "
            "surface because the extra ramification point collides with a residual "
            "pole; it is intentionally absent from the open input replay."
        ),
        "limitations": [
            "The residual quartic squarefree/coprime conditions are separate open conditions.",
            "Finite-field points need not lift to characteristic zero.",
            "No Mordell-Weil section is constructed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.surface)
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
        "finite_primes": [item["prime"] for item in payload["finite_field_replays"]],
        "finite_record_counts": [item["record_count"] for item in payload["finite_field_replays"]],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
