#!/usr/bin/env python3
"""Unique bounded CRT reconstruction on the exact open split IV loci.

The inputs are the exact finite-field loci with configuration

    I12 + I4 + IV + 4 I1

at the good primes 11, 17, and 19.  Every one of the 4*5*4 triples is matched
coordinatewise in

    p0,p1,p2,p3,r,s.

For M=11*17*19=3553 and B=42, 2*B^2=3528<M, so each residue class has at
most one reduced representation a/b with |a|<=B and 1<=b<=B.  A candidate is
promoted only after exact substitution over Q into all five surface equations,
exact fibre-order checks, both rational split-square conditions, and exact
residual-quartic squarefreeness/coprimality.

Absence within this box is not a global obstruction.  No section or rank-30
claim follows from a reconstructed surface alone.
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
    value = residues[0] % moduli[0]
    modulus = moduli[0]
    for residue, next_modulus in zip(residues[1:], moduli[1:], strict=True):
        if math.gcd(modulus, next_modulus) != 1:
            raise ValueError(("noncoprime CRT moduli", modulus, next_modulus))
        step = (
            (residue - value) * pow(modulus, -1, next_modulus)
        ) % next_modulus
        value = (value + modulus * step) % (modulus * next_modulus)
        modulus *= next_modulus
    return value


def fraction_mod_prime(value: Fraction, prime: int) -> int:
    denominator = value.denominator % prime
    if denominator == 0:
        raise ZeroDivisionError((value, prime))
    return (
        value.numerator % prime
        * pow(denominator, -1, prime)
    ) % prime


def load_open(path: Path, expected_prime: int) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if int(data["prime"]) != expected_prime:
        raise ValueError((path, data.get("prime"), expected_prime))
    expected_count = int(data["accepted_exact_configuration_count"])
    representatives = data["representatives"]
    if len(representatives) != expected_count:
        raise AssertionError(("representative count mismatch", path))
    for representative in representatives:
        parameters = representative["parameters"]
        if "e0" in parameters and int(parameters["e0"]) != 1:
            raise AssertionError(("noncanonical e0 representative", path))
    return data


def parameter_record(candidate, names):
    parameters = candidate["parameters"]
    return {name: int(parameters[name]) for name in names}


def build(
    *,
    f11_path: Path,
    f17_path: Path,
    f19_path: Path,
    bound: int,
) -> dict[str, object]:
    helper = load_helper()
    names = tuple(helper.PARAMETER_NAMES)
    primes = [11, 17, 19]
    loci = [
        load_open(f11_path, 11),
        load_open(f17_path, 17),
        load_open(f19_path, 19),
    ]
    modulus = math.prod(primes)
    if 2 * bound * bound >= modulus:
        raise ValueError((
            "bound does not ensure unique symmetric reconstruction",
            bound,
            modulus,
        ))

    exact_records: dict[str, dict[str, object]] = {}
    triple_records = []
    full_coordinate_candidate_count = 0

    ranges = [range(len(locus["representatives"])) for locus in loci]
    for indices in itertools.product(*ranges):
        candidates = [
            locus["representatives"][index]
            for locus, index in zip(loci, indices, strict=True)
        ]
        residues = [parameter_record(candidate, names) for candidate in candidates]
        crt_residues = {
            name: crt_many([record[name] for record in residues], primes)
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
                "unique reconstruction inequality was violated",
                indices,
                options,
            ))
        all_coordinates = all(len(options[name]) == 1 for name in names)
        parameters = None
        status = "NOT_ALL_COORDINATES_RECONSTRUCTED"
        promoted = []
        if all_coordinates:
            full_coordinate_candidate_count += 1
            parameters = {name: options[name][0] for name in names}
            for name in names:
                for prime, record in zip(primes, residues, strict=True):
                    if fraction_mod_prime(parameters[name], prime) != record[name] % prime:
                        raise AssertionError((
                            "reconstructed coordinate fails source reduction",
                            name,
                            parameters[name],
                            prime,
                            record[name],
                        ))
            exact = helper.verify_over_q(parameters)
            if exact is None:
                status = "EXACT_Q_SUBSTITUTION_REJECTED"
            else:
                status = "EXACT_Q_OPEN_SURFACE_VERIFIED"
                exact_records[exact["record_sha256"]] = exact
                promoted.append(exact["record_sha256"])

        triple_records.append({
            "representative_indices_f11_f17_f19": list(indices),
            "source_parameter_residues": {
                str(prime): record
                for prime, record in zip(primes, residues, strict=True)
            },
            "crt_residues_mod_3553": crt_residues,
            "bounded_rational_options": {
                name: [str(value) for value in options[name]]
                for name in names
            },
            "all_coordinates_reconstruct": all_coordinates,
            "reconstructed_parameters": (
                None if parameters is None else
                {name: str(parameters[name]) for name in names}
            ),
            "exact_verification_status": status,
            "promoted_record_sha256": promoted,
        })

    exact = [exact_records[key] for key in sorted(exact_records)]
    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_open_surface_crt_f11_f17_f19_unique_bounded",
        "truth_status": (
            "EXACT unique bounded rational reconstruction and Q substitution "
            "from the open split F11/F17/F19 loci; absence outside the declared "
            "box is unresolved; no section or rank-30 conclusion"
        ),
        "source_open_locus_record_sha256": {
            str(prime): locus["record_sha256"]
            for prime, locus in zip(primes, loci, strict=True)
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
        "open_surface_counts": {
            str(prime): locus["accepted_exact_configuration_count"]
            for prime, locus in zip(primes, loci, strict=True)
        },
        "surface_triple_count": len(triple_records),
        "full_coordinate_candidate_count": full_coordinate_candidate_count,
        "exact_q_open_surface_count": len(exact),
        "exact_q_open_surfaces": exact,
        "triples": triple_records,
        "limitations": [
            "Coordinatewise CRT does not identify a characteristic-zero point unless exact Q substitution succeeds.",
            "Only reduced rationals with bounded numerator and denominator are tested.",
            "A Q-point with bad or nonintegral reduction at one source prime is invisible to this gate.",
            "Other projective charts or equivalent normalizations are not covered.",
            "A rational IV surface without the required height-79/12 section is insufficient.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f11", type=Path, required=True)
    parser.add_argument("--f17", type=Path, required=True)
    parser.add_argument("--f19", type=Path, required=True)
    parser.add_argument("--bound", type=int, default=42)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(
        f11_path=arguments.f11,
        f17_path=arguments.f17,
        f19_path=arguments.f19,
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
        "exact_q_open_surface_count": payload["exact_q_open_surface_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
