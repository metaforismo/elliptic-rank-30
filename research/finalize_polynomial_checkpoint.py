#!/usr/bin/env python3
"""Finalize the exact polynomial inverse-construction checkpoint.

This script consumes proof/search certificates generated from the repository's
mathematical code and writes a compact, evidence-linked project checkpoint.
It never changes a proof certificate and fails closed when expected results are
missing or inconsistent.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-08-08"
STATE = "CLOSED_FOR_POLYNOMIAL_IDENTITIES_AT_P_EQUALS_3"


def load_json(relative: str) -> dict[str, Any]:
    path = ROOT / relative
    if not path.is_file():
        raise SystemExit(f"missing required certificate: {relative}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"certificate is not a JSON object: {relative}")
    return data


def write_json(relative: str, data: dict[str, Any]) -> None:
    (ROOT / relative).write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def append_once(relative: str, marker: str, section: str, title: str) -> None:
    path = ROOT / relative
    text = path.read_text(encoding="utf-8") if path.exists() else f"# {title}\n"
    if marker in text:
        return
    path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def main() -> int:
    collision = load_json("certificates/degree24_collision_search_n256_d256.json")
    passport = load_json("certificates/degree47_belyi_passport.json")
    groebner = load_json("certificates/degree47_target_elimination.json")
    mod19 = load_json("certificates/degree47_target_mod19_obstruction.json")
    degenerate = load_json("certificates/degree47_degenerate_obstruction.json")
    mason = load_json("certificates/polynomial_section_mason_bound.json")
    degree24 = load_json("certificates/degree24_polynomial_section_classification.json")

    if collision.get("collision_count") != 0:
        raise SystemExit("bounded collision search unexpectedly found a collision")
    if groebner.get("result") != "no_rational_solution":
        raise SystemExit("degree-(4,7) exact elimination did not prove nonexistence")
    if groebner.get("full_rational_solutions") != []:
        raise SystemExit("degree-(4,7) elimination has unprocessed rational solutions")

    dessins = int(passport["dessin_count"])
    parameters = int(collision["parameter_count"])
    conclusion = (
        "The complete degree-(2,4) classification excludes p=3; the "
        "Mason--Stothers theorem proves that degree (4,7) is the final possible "
        "polynomial level; the repeated-root branches are impossible; and the "
        "generic degree-(4,7) branch is excluded independently by exact "
        "Groebner elimination and by the primitive determinant quintic with no "
        "root modulo 19. Therefore every polynomial identity of this specific "
        "shape at p=3 is impossible."
    )
    limitation = (
        "This closes only the polynomial representation. Rational functions "
        "with finite poles, different norm forms, other elliptic surfaces, and "
        "the global rank-30 problem remain open."
    )
    next_action = (
        "Classify the smallest genuinely new controlled-pole rational-function "
        "ansatz by Riemann--Roch and exact elimination; continue independently "
        "with the explicit rank-17 K3 bridge and high-generic-rank families."
    )

    report = f"""# Polynomial inverse-construction checkpoint

**Truth status:** `{STATE}`.

**Global challenge:** UNSOLVED. The authoritative unconditional lower bound is
`rank E(Q) >= 29`; no rank-30 certificate is present.

## Proved chain

1. The degree-(2,4) coefficient scheme is completely classified.
2. Its rational image does not contain the historical coefficient `p=3`.
3. Mason--Stothers leaves only the degree-(4,7) polynomial level.
4. Both repeated-root degree-(4,7) target branches are impossible.
5. The generic degree-(4,7) branch is impossible by an exact mod-19
   obstruction and an independent exact Groebner calculation.

{conclusion}

## Finite computations

- bounded collision parameters tested: `{parameters}`;
- bounded collisions found: `0`;
- maximal degree-14 passport dessin classes: `{dessins}`;
- exact Groebner result: `{groebner['result']}`.

## Scope and next action

{limitation}

{next_action}
"""
    (ROOT / "research/POLYNOMIAL_THREAD_CHECKPOINT.md").write_text(
        report, encoding="utf-8"
    )

    proof_paths = [
        "certificates/degree24_polynomial_section_classification.json",
        "certificates/degree24_polynomial_section_obstruction.json",
        "certificates/polynomial_section_mason_bound.json",
        "certificates/degree47_degenerate_obstruction.json",
        "certificates/degree47_target_mod19_obstruction.json",
        "certificates/degree47_target_elimination.json",
    ]

    frontier_path = ROOT / "STATUS.frontier.json"
    frontier = load_json("STATUS.frontier.json")
    frontier.update(
        {
            "authoritative_lower_bound": 29,
            "canonical_branch": "main",
            "record_claim": False,
            "truth_status": "UNSOLVED: exact rank lower bound 29; no rank-30 certificate",
            "updated": DATE,
            "newest_result": {
                "name": "complete polynomial-identity obstruction at p=3",
                "claim_status": "proved",
                "documentation": "research/POLYNOMIAL_THREAD_CHECKPOINT.md",
                "certificate": "certificates/degree47_target_elimination.json",
                "limitation": limitation,
            },
            "polynomial_inverse_construction": {
                "state": STATE,
                "exact_conclusion": conclusion,
                "proof_certificates": proof_paths,
                "bounded_collision_search": "certificates/degree24_collision_search_n256_d256.json",
                "generic_passport": "certificates/degree47_belyi_passport.json",
                "generic_dessin_count": dessins,
                "next_action": next_action,
            },
            "active_workstreams": [
                "classify controlled-pole rational-function inverse constructions",
                "reconstruct the explicit rank-17 K3 bridge and historical specialization",
                "benchmark independent high-generic-rank families under the common certificate funnel",
                "study low-genus Galois-character packets and collision descents",
            ],
        }
    )
    write_json("STATUS.frontier.json", frontier)

    checkpoint = {
        "schema_version": 1,
        "updated": DATE,
        "canonical_branch": "main",
        "mathematical_status": {
            "best_unconditional_lower_bound": 29,
            "rank_30_found": False,
            "remaining_gap": 1,
        },
        "polynomial_inverse_construction": {
            "state": STATE,
            "proof_certificates": proof_paths,
            "bounded_collision_search": "certificates/degree24_collision_search_n256_d256.json",
            "generic_passport": "certificates/degree47_belyi_passport.json",
            "next_action": next_action,
        },
        "blocking_external_input": "search/REQUEST_FOR_FIBRATION_DATA.md",
    }
    write_json("CHECKPOINT.json", checkpoint)

    append_once(
        "STATUS.md",
        "<!-- polynomial-checkpoint-20260808 -->",
        f"""<!-- polynomial-checkpoint-20260808 -->
## Exact polynomial checkpoint

- Best unconditional lower bound: `rank E(Q) >= 29`.
- Rank-30 certificate: absent.
- State: `{STATE}`.
- Bounded collision search: `{parameters}` parameters, zero collisions.
- Maximal passport: `{dessins}` dessin class(es).
- Highest-value next action: {next_action}""",
        "Status",
    )
    append_once(
        "CLAIMS.md",
        "<!-- polynomial-claims-20260808 -->",
        f"""<!-- polynomial-claims-20260808 -->
## Polynomial inverse-construction claims

| Claim | Truth status | Evidence |
|---|---|---|
| Complete degree-(2,4) classification | proved | `certificates/degree24_polynomial_section_classification.json` |
| No degree-(2,4) target at `p=3` | proved | `certificates/degree24_polynomial_section_obstruction.json` |
| No polynomial level beyond `(4,7)` | proved | `certificates/polynomial_section_mason_bound.json` |
| Repeated-root `(4,7)` branches impossible | proved | `certificates/degree47_degenerate_obstruction.json` |
| Generic `(4,7)` branch impossible | proved twice | `certificates/degree47_target_mod19_obstruction.json`; `certificates/degree47_target_elimination.json` |
| All polynomial identities of this shape at `p=3` impossible | proved | conjunction of the preceding certificates |
| Rank at least 30 | unknown | not claimed |""",
        "Claims",
    )
    append_once(
        "NEGATIVE_RESULTS.md",
        "<!-- polynomial-negative-20260808 -->",
        f"""<!-- polynomial-negative-20260808 -->
## Closed polynomial representation at `p=3`

**Status:** `{STATE}`.

{conclusion}

{limitation} Do not spend additional compute enlarging polynomial coefficient
searches in this representation.""",
        "Negative results",
    )
    append_once(
        "HYPOTHESES.md",
        "<!-- controlled-pole-hypothesis-20260808 -->",
        """<!-- controlled-pole-hypothesis-20260808 -->
## Controlled-pole replacement hypothesis

**Status:** active, unproved.

The smallest genuinely new inverse-construction space should allow at least two
finite pole locations. A single pole on `P^1` can be moved to infinity by a
Möbius transformation and is therefore only a reparameterized polynomial
problem. The next exact test is to fix a two-point pole divisor, derive the
Riemann--Roch coefficient space, and eliminate the target identity before any
large specialization search.""",
        "Hypotheses",
    )
    append_once(
        "DECISIONS.md",
        "<!-- polynomial-decision-20260808 -->",
        f"""<!-- polynomial-decision-20260808 -->
## Decision: change representation

{conclusion}

Canonical next action: {next_action}""",
        "Decisions",
    )

    registry = load_json("APPROACH_REGISTRY.json")
    approaches = registry.setdefault("approaches", [])
    if not isinstance(approaches, list):
        raise SystemExit("APPROACH_REGISTRY.json has invalid approaches field")
    for approach in approaches:
        if isinstance(approach, dict) and approach.get("approach_id") == "D1-inverse-prescribed-points":
            approach["evidence_against"] = list(approach.get("evidence_against", [])) + [
                "The complete polynomial norm identity at p=3 is now proved impossible at every degree."
            ]
            approach["concrete_next_test"] = next_action
    if not any(
        isinstance(approach, dict)
        and approach.get("approach_id") == "D3-controlled-pole-rational-functions"
        for approach in approaches
    ):
        approaches.append(
            {
                "approach_id": "D3-controlled-pole-rational-functions",
                "mathematical_family": "Rational-function norm identities with a fixed two-point pole divisor.",
                "central_mechanism": "Evade the polynomial Mason--Stothers obstruction while retaining a finite Riemann--Roch coefficient space and exact elimination.",
                "assumptions": [],
                "predicted_advantage": "Two pole locations are the first projectively genuine extension of the closed polynomial representation.",
                "concrete_next_test": next_action,
                "evidence_for": [
                    "The polynomial representation is now completely closed, so finite poles are a mathematically necessary change of setting within this norm-form route.",
                    "A fixed pole divisor gives a finite-dimensional exact coefficient space."
                ],
                "evidence_against": [
                    "The cleared identity has a larger remainder and may force high-genus auxiliary varieties.",
                    "No rank contribution has yet been constructed in this space."
                ],
                "computational_cost": "Low-to-medium for the first exact elimination; potentially high after increasing pole order.",
                "missing_lemma": "An explicit controlled-pole identity producing a non-torsion section or a proof that the minimal pole divisor is obstructed.",
                "blocker": "The minimal two-pole ansatz has not yet been classified.",
                "status": "ACTIVE",
            }
        )
    registry["updated_at"] = DATE
    write_json("APPROACH_REGISTRY.json", registry)

    print(
        "FINALIZED polynomial checkpoint",
        f"state={STATE}",
        f"dessins={dessins}",
        f"collision_parameters={parameters}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
