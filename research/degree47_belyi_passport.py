#!/usr/bin/env python3
"""Enumerate transitive dessins for the maximal degree-(4,7) passport.

The passport is

    sigma_0:      2^7
    sigma_1:      3^4 2
    sigma_infty:  11 1^3,

with sigma_0 sigma_1 sigma_infty = 1.  We fix sigma_infty to a canonical
11-cycle plus three fixed points, enumerate every perfect matching sigma_0,
and quotient by the full centralizer C_11 x S_3 of sigma_infty.

The output is an exact finite combinatorial certificate.  It does not by
itself prove that the corresponding Belyi maps are defined over Q.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Tuple

Permutation = Tuple[int, ...]


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def compose(left: Permutation, right: Permutation) -> Permutation:
    if len(left) != len(right):
        raise ValueError("permutations have different degrees")
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def conjugate(permutation: Permutation, by: Permutation) -> Permutation:
    return compose(by, compose(permutation, inverse(by)))


def cycle_partition(permutation: Permutation) -> Tuple[int, ...]:
    seen = [False] * len(permutation)
    cycles: list[int] = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        current = start
        length = 0
        while not seen[current]:
            seen[current] = True
            length += 1
            current = permutation[current]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def is_transitive(generators: Sequence[Permutation]) -> bool:
    degree = len(generators[0])
    reached = {0}
    frontier = [0]
    inverses = [inverse(generator) for generator in generators]
    while frontier:
        point = frontier.pop()
        for generator in itertools.chain(generators, inverses):
            image = generator[point]
            if image not in reached:
                reached.add(image)
                frontier.append(image)
    return len(reached) == degree


def canonical_sigma_infinity() -> Permutation:
    return tuple(list(range(1, 11)) + [0, 11, 12, 13])


def perfect_matchings(points: Tuple[int, ...]) -> Iterator[Permutation]:
    degree = len(points)
    if degree % 2:
        raise ValueError("a perfect matching needs even degree")
    partner = [-1] * degree

    def recurse(remaining: Tuple[int, ...]) -> Iterator[Permutation]:
        if not remaining:
            yield tuple(partner)
            return
        first = remaining[0]
        for offset in range(1, len(remaining)):
            second = remaining[offset]
            partner[first] = second
            partner[second] = first
            next_remaining = remaining[1:offset] + remaining[offset + 1 :]
            yield from recurse(next_remaining)
            partner[first] = -1
            partner[second] = -1

    yield from recurse(points)


def centralizer_of_sigma_infinity() -> Tuple[Permutation, ...]:
    degree = 14
    result: list[Permutation] = []
    for rotation in range(11):
        for fixed_permutation in itertools.permutations((11, 12, 13)):
            mapping = [0] * degree
            for point in range(11):
                mapping[point] = (point + rotation) % 11
            for source, target in zip((11, 12, 13), fixed_permutation):
                mapping[source] = target
            result.append(tuple(mapping))
    if len(result) != 66 or len(set(result)) != 66:
        raise AssertionError("centralizer enumeration failed")
    sigma_infinity = canonical_sigma_infinity()
    if any(compose(element, sigma_infinity) != compose(sigma_infinity, element) for element in result):
        raise AssertionError("a purported centralizer element does not commute")
    return tuple(result)


def canonical_representative(
    sigma_0: Permutation,
    centralizer: Sequence[Permutation],
) -> Permutation:
    return min(conjugate(sigma_0, element) for element in centralizer)


def permutation_cycles(permutation: Permutation) -> list[list[int]]:
    seen = [False] * len(permutation)
    result: list[list[int]] = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cycle: list[int] = []
        current = start
        while not seen[current]:
            seen[current] = True
            cycle.append(current)
            current = permutation[current]
        result.append(cycle)
    return sorted(result, key=lambda cycle: (-len(cycle), cycle))


def enumerate_dessins() -> dict:
    degree = 14
    sigma_infinity = canonical_sigma_infinity()
    sigma_infinity_inverse = inverse(sigma_infinity)
    centralizer = centralizer_of_sigma_infinity()
    target_sigma_1 = (3, 3, 3, 3, 2)

    raw_count = 0
    transitive_count = 0
    representatives: dict[Permutation, dict] = {}

    for sigma_0 in perfect_matchings(tuple(range(degree))):
        sigma_1 = compose(sigma_0, sigma_infinity_inverse)
        if cycle_partition(sigma_1) != target_sigma_1:
            continue
        raw_count += 1
        if not is_transitive((sigma_0, sigma_1)):
            continue
        transitive_count += 1
        canonical = canonical_representative(sigma_0, centralizer)
        if canonical in representatives:
            continue
        canonical_sigma_1 = compose(canonical, sigma_infinity_inverse)
        automorphisms = sum(
            conjugate(canonical, element) == canonical for element in centralizer
        )
        triple_check = compose(canonical, compose(canonical_sigma_1, sigma_infinity))
        if triple_check != identity(degree):
            raise AssertionError("permutation triple relation failed")
        representatives[canonical] = {
            "automorphism_group_order": automorphisms,
            "sigma_0": permutation_cycles(canonical),
            "sigma_1": permutation_cycles(canonical_sigma_1),
            "sigma_infinity": permutation_cycles(sigma_infinity),
        }

    records = [representatives[key] for key in sorted(representatives)]
    result = {
        "certificate_id": "degree47-belyi-passport-enumeration-v1",
        "centralizer_order": len(centralizer),
        "claim_status": "proved finite combinatorial enumeration",
        "degree": degree,
        "dessin_count": len(records),
        "dessins": records,
        "fixed_sigma_infinity": permutation_cycles(sigma_infinity),
        "passport": {
            "0": [2, 2, 2, 2, 2, 2, 2],
            "1": [3, 3, 3, 3, 2],
            "infinity": [11, 1, 1, 1],
        },
        "raw_matching_count": 135135,
        "raw_passport_count": raw_count,
        "schema_version": 1,
        "scope_limitation": "Enumeration of topological dessins only; fields of definition and exact Belyi maps require subsequent arithmetic reconstruction.",
        "transitive_passport_count": transitive_count,
        "verification": {
            "equivalence": "simultaneous conjugacy by the full centralizer C_11 x S_3 of the fixed sigma_infinity",
            "method": "exhaustive enumeration of all 13!! perfect matchings",
            "relation": "sigma_0 sigma_1 sigma_infinity = identity",
        },
    }
    unhashed = json.dumps(result, sort_keys=True, separators=(",", ":")).encode("utf-8")
    result["record_sha256"] = hashlib.sha256(unhashed).hexdigest()
    return result


def canonical_payload(result: dict) -> str:
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()

    result = enumerate_dessins()
    payload = canonical_payload(result)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    if args.check and args.check.read_text(encoding="utf-8") != payload:
        raise SystemExit(f"certificate mismatch: {args.check}")
    print(
        "VERIFIED maximal Belyi passport enumeration",
        f"dessins={result['dessin_count']}",
        f"transitive_labelled={result['transitive_passport_count']}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
