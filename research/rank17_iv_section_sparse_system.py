#!/usr/bin/env sage -python
"""Build a sparse exact section-incidence system for one additive-IV surface.

The endpoint recurrence used by earlier probes removes z_1,...,z_9, but the
resulting rational numerators have very high degree.  This builder keeps those
nine square-root coefficients explicit.  For fixed q != 0 it produces an
exact finite-field system with

    d0, d1, w1,...,w6, z1,...,z9, inv

and the coefficient identities for Z^2=H, the selected IV tangent, and a
Rabinowitsch equation excluding d0*D(1)*z0=0.  It intentionally postpones the
resultant(D,X) filter; every returned point must be filtered for minimality.

The output is suitable for msolve's documented text interface.  This script
makes no p-adic, characteristic-zero, or rank-30 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.version import version as sage_version


def convolution(left, right):
    zero = left[0] * 0
    output = [zero for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i + j] += a * b
    return output


def divide_by_t_minus_one_squared(coefficients, ring):
    remainder = list(coefficients)
    quotient = [ring(0) for _ in range(len(coefficients) - 2)]
    for degree in range(len(coefficients) - 1, 1, -1):
        coefficient = remainder[degree]
        quotient[degree - 2] = coefficient
        remainder[degree - 2] -= coefficient
        remainder[degree - 1] += 2 * coefficient
        remainder[degree] -= coefficient
    if remainder[0] != 0 or remainder[1] != 0:
        raise AssertionError("invariant is not divisible by (t-1)^2")
    return quotient


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_candidate(path: Path, index: int):
    certificate = json.loads(path.read_text(encoding="utf-8"))
    candidates = certificate["candidates"]
    if index < 0 or index >= len(candidates):
        raise IndexError((index, len(candidates)))
    return certificate, candidates[index]


def build_system(
    *,
    certificate_path: Path,
    candidate_index: int,
    q_value: int,
    i4_root: int,
    iv_root: int,
):
    certificate, candidate = load_candidate(certificate_path, candidate_index)
    prime = int(certificate["prime"])
    field = GF(prime)
    q = field(q_value)
    if q == 0:
        raise ValueError("this sparse chart requires q != 0")
    if q**2 + 3 == 0:
        raise ValueError("q parametrizes the singular infinity point")

    parameters = candidate["parameters"]
    e0 = field(parameters["e0"])
    p0 = field(parameters["p0"])
    p1 = field(parameters["p1"])
    rho = field(i4_root)
    eta = field(iv_root)

    valid_i4 = {
        int(value)
        for value in candidate["split_tangent_checks"]["i4_square_roots"]
    }
    valid_iv = {
        int(value)
        for value in candidate["split_tangent_checks"]["iv_square_roots"]
    }
    if int(rho) not in valid_i4:
        raise ValueError(("invalid I4 tangent root", int(rho), sorted(valid_i4)))
    if int(eta) not in valid_iv:
        raise ValueError(("invalid IV tangent root", int(eta), sorted(valid_iv)))
    if rho**2 != -3 * e0 * p0:
        raise AssertionError("I4 tangent square check failed")
    b_at_one = field(
        candidate["exact_fibre_checks"]["c6_quadratic_coefficient_at_one"]
    )
    if eta**2 != -2 * b_at_one:
        raise AssertionError("IV tangent square check failed")

    names = (
        "d0", "d1", "w1", "w2", "w3", "w4", "w5", "w6",
        "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9",
        "inv",
    )
    ring = PolynomialRing(field, names=names, order="degrevlex")
    (
        d0, d1, w1, w2, w3, w4, w5, w6,
        z1, z2, z3, z4, z5, z6, z7, z8, z9,
        inv,
    ) = ring.gens()

    c4 = [ring(field(value)) for value in candidate["c4_coefficients_ascending"]]
    c6 = [ring(field(value)) for value in candidate["c6_coefficients_ascending"]]
    A = divide_by_t_minus_one_squared(c4, ring)
    B = divide_by_t_minus_one_squared(c6, ring)

    D = [d0, d1, ring(1)]
    D2 = convolution(D, D)
    D4 = convolution(D2, D2)
    D6 = convolution(D4, D2)

    w0 = e0 * p0 * d0**2
    w7 = ring(q**2 + 2)
    W = [w0, w1, w2, w3, w4, w5, w6, w7]
    W2 = convolution(W, W)
    W3 = convolution(W2, W)
    A_W_D4 = convolution(convolution(A, W), D4)
    core = [W3[index] - 3 * A_W_D4[index] for index in range(22)]
    times_t_minus_one = [ring(0) for _ in range(23)]
    for index, coefficient in enumerate(core):
        times_t_minus_one[index] -= coefficient
        times_t_minus_one[index + 1] += coefficient
    B_D6 = convolution(B, D6)
    bracket = [
        times_t_minus_one[index] - 2 * B_D6[index]
        for index in range(23)
    ]
    if bracket[0] != 0 or bracket[1] != 0:
        raise AssertionError("I4-adapted numerator is not divisible by t^2")
    H = bracket[2:]
    if len(H) != 21:
        raise AssertionError(("unexpected square-target length", len(H)))

    x1 = w0 - w1
    moving_node_first_jet = x1 + e0 * (p0 * D2[1] + p1 * D2[0])
    z0 = -rho * d0 * moving_node_first_jet
    z10 = ring(q * (q**2 + 3))
    z = [z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10]

    if H[0] - z0**2 != 0:
        raise AssertionError("constant square identity is not automatic")
    if H[20] - z10**2 != 0:
        raise AssertionError("leading square identity is not automatic")

    equations = []
    labels = []
    for degree in range(1, 20):
        square_coefficient = ring(0)
        for index in range(11):
            other = degree - index
            if 0 <= other < 11:
                square_coefficient += z[index] * z[other]
        equations.append(square_coefficient - H[degree])
        labels.append(f"square_coefficient_{degree}")

    D1 = d0 + d1 + 1
    equations.append(sum(z, ring(0)) - eta * D1**3)
    labels.append("chosen_IV_tangent")

    saturation = d0 * D1 * z0
    equations.append(inv * saturation - 1)
    labels.append("rabinowitsch_d0_D1_z0")

    metadata = {
        "schema_version": 1,
        "truth_status": (
            f"Exact sparse section-incidence system over F_{prime}; "
            "the resultant(D,X) minimality filter is postponed; no p-adic, "
            "characteristic-zero, or rank-30 conclusion."
        ),
        "sage_version": str(sage_version),
        "prime": prime,
        "surface_certificate_sha256": certificate.get("certificate_sha256"),
        "candidate_index": candidate_index,
        "candidate_record_sha256": candidate.get("record_sha256"),
        "q": int(q),
        "i4_root": int(rho),
        "iv_root": int(eta),
        "variables": list(names),
        "equation_labels": labels,
        "equation_count": len(equations),
        "equation_total_degrees": [int(value.total_degree()) for value in equations],
        "equation_term_counts": [len(value.dict()) for value in equations],
        "automatic_square_degrees": [0, 20],
        "postponed_filters": ["resultant_t(D,X) != 0"],
    }
    metadata["system_sha256"] = canonical_hash({
        "variables": list(names),
        "prime": prime,
        "equations": [str(value) for value in equations],
    })
    return ring, equations, metadata


def write_msolve(path: Path, ring, equations, prime: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(str(value) for value in ring.gens()), str(prime)]
    for index, equation in enumerate(equations):
        suffix = "," if index + 1 < len(equations) else ""
        lines.append(f"{equation}{suffix}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--candidate-index", type=int, required=True)
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--i4-root", type=int, required=True)
    parser.add_argument("--iv-root", type=int, required=True)
    parser.add_argument("--msolve-input", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    arguments = parser.parse_args()

    ring, equations, metadata = build_system(
        certificate_path=arguments.certificate,
        candidate_index=arguments.candidate_index,
        q_value=arguments.q,
        i4_root=arguments.i4_root,
        iv_root=arguments.iv_root,
    )
    write_msolve(arguments.msolve_input, ring, equations, metadata["prime"])
    metadata["msolve_input_sha256"] = hashlib.sha256(
        arguments.msolve_input.read_bytes()
    ).hexdigest()
    arguments.metadata.parent.mkdir(parents=True, exist_ok=True)
    arguments.metadata.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
