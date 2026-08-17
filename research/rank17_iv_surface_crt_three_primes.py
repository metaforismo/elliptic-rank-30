#!/usr/bin/env python3
"""Unique bounded rational reconstruction from the F11/F13/F17 IV loci.

The input files are exact geometric-quotient certificates.  Every triple of
finite-field surfaces is matched coordinatewise in the six normalized
parameters

    p0, p1, p2, p3, r, s.

For M=11*13*17=2431 and B=34, the inequality

    2*B^2 < M

ensures that each residue has at most one reduced representation a/b with
|a|<=B and 1<=b<=B.  A six-coordinate candidate is promoted only after exact
substitution into all five surface equations, exact discriminant-order checks,
the two rational split-square conditions, and squarefreeness/coprimality of
the residual quartic.

Failure inside this declared box is not a global obstruction and says nothing
about rational sections or rank 30.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER_PATH = HERE / "rank17_iv_surface_crt_reconstruction.py"


def load_helper():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_surface_crt_reconstruction", HELPER_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(HELPER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def crt_many(residues: list[int], moduli: list[int]) -> int:
    if len(residues) != len(moduli) or not residues:
        raise ValueError("CRT input dimensions differ or are empty")
    value = residues[0] % moduli[0]
    modulus = moduli[0]
    for residue, next_modulus in zip(residues[1:], moduli[1:], strict=True):
        if math.gcd(modulus, next_modulus) != 1:
            raise ValueError(("noncoprime CRT moduli", modulus, next_modulus))
        step = (
            (residue - value) * pow(modulus, -1, next_modulus)
        ) % next_modulus
        value += modulus * step
        modulus *= next_modulus
        value %= modulus
    return value


def fraction_mod_prime(value: Fraction, prime: int) -> int:
    denominator = value.denominator % prime
    if denominator == 0:
        raise ZeroDivisionError((value, prime))
    return (
        value.numerator % prime
        * pow(denominator, -1, prime)
    ) % prime


def load_quotient(path: Path, expected_prime: int) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if int(data["prime"]) != expected_prime:
        raise ValueError((path, data.get("prime"), expected_prime))
    if 2 * int(data["geometric_surface_count"]) != int(data["raw_candidate_count"]):
        raise AssertionError(("invalid sign quotient", path))
    if len(data["representatives"]) != int(data["geometric_surface_count"]):
        raise AssertionError(("representative count mismatch", path))
    for representative in data["representatives"]:
        if int(representative["parameters"]["e0"]) != 1:
            raise AssertionError(("noncanonical e0 representative", path))
    return data


def parameter_residues(candidate: dict[str, object], names) -> dict[str, int]:
    parameters = candidate["parameters"]
    return {name: int(parameters[name]) for name in names}


def build(
    *,
    f11_path: Path,
    f13_path: Path,
    f17_path: Path,
    bound: int,
) -> dict[str, object]:
    helper = load_helper()
    names = tuple(helper.PARAMETER_NAMES)
    primes = [11, 13, 17]
    quotient11 = load_quotient(f11_path, 11)
    quotient13 = load_quotient(f13_path, 13)
    quotient17 = load_quotient(f17_path, 17)
    quotients = [quotient11, quotient13, quotient17]
    modulus = math.prod(primes)
    if 2 * bound * bound >= modulus:
        raise ValueError((
            "bound does not give unique symmetric rational reconstruction",
            bound,
            modulus,
        ))

    exact_records: dict[str, dict[str, object]] = {}
    triple_records = []
    full_coordinate_candidate_count = 0

    products = itertools.product(
        range(len(quotient11["representatives"])),
        range(len(quotient13["representatives"])),
        range(len(quotient17["representatives"])),
    )
    for indices in products:
        candidates = [
            quotient["representatives"][index]
            for quotient, index in zip(quotients, indices, strict=True)
        ]
        residue_records = [
            parameter_residues(candidate, names)
            for candidate in candidates
        ]
        crt_residues = {
            name: crt_many(
                [record[name] for record in residue_records],
                primes,
            )
            for name in names
        }
        options = {
            name: helper.bounded_rational_options(
                crt_residues[name], modulus, bound
            )
            for name in names
        }
        if any(len(values) > 1 for values in options.values()):
            raise AssertionError((
                "rational reconstruction uniqueness theorem was violated",
                indices,
                options,
            ))
        all_coordinates_reconstruct = all(
            len(options[name]) == 1 for name in names
        )
        promoted = []
        parameters = None
        verification_status = "NOT_ALL_COORDINATES_RECONSTRUCTED"
        if all_coordinates_reconstruct:
            full_coordinate_candidate_count += 1
            parameters = {
                name: options[name][0]
                for name in names
            }
            for name in names:
                for prime, residue in zip(primes, residue_records, strict=True):
                    if fraction_mod_prime(parameters[name], prime) != residue[name] % prime:
                        raise AssertionError((
                            "reconstructed coordinate fails reduction",
                            name,
                            parameters[name],
                            prime,
                            residue[name],
                        ))
            exact = helper.verify_over_q(parameters)
            if exact is None:
                verification_status = "EXACT_Q_SUBSTITUTION_REJECTED"
            else:
                verification_status = "EXACT_Q_SURFACE_VERIFIED"
                exact_records[exact["record_sha256"]] = exact
                promoted.append(exact["record_sha256"])

        triple_records.append({
            "representative_indices_f11_f13_f17": list(indices),
            "source_parameter_residues": {
                str(prime): residue
                for prime, residue in zip(primes, residue_records, strict=True)
            },
            "crt_residues_mod_2431": crt_residues,
            "bounded_rational_options": {
                name: [str(value) for value in options[name]]
                for name in names
            },
            "all_coordinates_reconstruct": all_coordinates_reconstruct,
            "reconstructed_parameters": (
                None if parameters is None else
                {name: str(parameters[name]) for name in names}
            ),
            "exact_verification_status": verification_status,
            "promoted_record_sha256": promoted,
        })

    exact = [exact_records[key] for key in sorted(exact_records)]
    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_surface_crt_f11_f13_f17_unique_bounded",
        "truth_status": (
            "EXACT unique bounded rational reconstruction and Q substitution "
            "from the geometric F11/F13/F17 surface loci; absence outside the "
            "declared bound is unresolved; no section or rank-30 conclusion"
        ),
        "source_quotient_record_sha256": {
            "11": quotient11["record_sha256"],
            "13": quotient13["record_sha256"],
            "17": quotient17["record_sha256"],
        },
        "primes": primes,
        "modulus": modulus,
        "parameter_names": list(names),
        "rational_numerator_denominator_bound": bound,
        "uniqueness_criterion": {
            "formula": "2*B^2 < M",
            "left_hand_side": 2 * bound * bound,
            "right_hand_side": modulus,
            "satisfied": True,
        },
        "geometric_surface_counts": {
            "11": quotient11["geometric_surface_count"],
            "13": quotient13["geometric_surface_count"],
            "17": quotient17["geometric_surface_count"],
        },
        "surface_triple_count": len(triple_records),
        "full_coordinate_candidate_count": full_coordinate_candidate_count,
        "exact_q_surface_count": len(exact),
        "exact_q_surfaces": exact,
        "triples": triple_records,
        "limitations": [
            "Coordinatewise CRT does not prove that unrelated finite-field points belong to one characteristic-zero point.",
            "Only reduced rationals with |numerator| and denominator at most the declared bound are tested.",
            "The affine normalization may omit other projective charts or equivalent models.",
            "A rational IV surface without the required height-79/12 section is insufficient.",
            "This certificate does not change the unconditional rank-29 lower bound."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--f11", type=Path, required=True)
    parser.add_argument("--f13", type=Path, required=True)
    parser.add_argument("--f17", type=Path, required=True)
    parser.add_argument("--bound", type=int, default=34)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()

    payload = build(
        f11_path=arguments.f11,
        f13_path=arguments.f13,
        f17_path=arguments.f17,
        bound=arguments.bound,
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
        "surface_triple_count": payload["surface_triple_count"],
        "full_coordinate_candidate_count": payload["full_coordinate_candidate_count"],
        "exact_q_surface_count": payload["exact_q_surface_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
