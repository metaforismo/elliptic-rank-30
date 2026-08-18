#!/usr/bin/env python3
"""Combine the three exact certificates that close the split additive-IV branch."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_and_verify(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = data.get("record_sha256")
    if not expected:
        raise AssertionError(("source certificate has no record hash", path))
    body = {key: value for key, value in data.items() if key != "record_sha256"}
    actual = canonical_hash(body)
    if actual != expected:
        raise AssertionError(("source record hash mismatch", path, expected, actual))
    return data


def build(singular_path: Path, open_path: Path, split_path: Path) -> dict[str, object]:
    singular = load_and_verify(singular_path)
    open_component = load_and_verify(open_path)
    split = load_and_verify(split_path)

    assert singular["solved_a"] is False
    assert singular["solved_b"] is False
    assert singular["unconditional_global_rank_lower_bound"] == 29
    assert singular["generic_q8_branch"]["resultant_rational_roots"] == {}
    assert singular["exceptional_q8_slope_zero"]["resultant_gcd"] == "1"

    assert open_component["solved_a"] is False
    assert open_component["solved_b"] is False
    assert open_component["unconditional_global_rank_lower_bound"] == 29
    assert open_component["residuals_contained_in_J"] is True
    assert open_component["component"]["plane_resultant_irreducible_over_QQ"] is True
    assert open_component["radical_membership"]["R1"]["rabinowitsch_unit_ideal"] is True
    assert open_component["radical_membership"]["R2"]["rabinowitsch_unit_ideal"] is True
    assert open_component["set_theoretic_consequence"] == (
        "V(I_open)=V(R1,R2)_open over the algebraic closure of Q."
    )

    assert split["solved_a"] is False
    assert split["solved_b"] is False
    assert split["unconditional_global_rank_lower_bound"] == 29
    assert split["plane_component"]["irreducible_over_Q"] is True
    assert split["plane_component"]["resultant_z_equals_quintic"] is True
    assert split["logarithmic_derivative_reconstruction"][
        "differential_equation_zero_over_Qw"
    ] is True
    assert split["split_square_classes"]["positive_set_intersection"] == "EmptySet"

    payload: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_split_branch_closure",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "classification": "ADDITIVE_IV_SPLIT_REALIZATION_CLOSED_OVER_Q",
        "truth_status": (
            "CERTIFIED over Q within the declared normalized logarithmic-derivative "
            "model: no nondegenerate surface with fibres I12+I4+IV has both the "
            "I4 and IV split conditions required by the additive rank-17 seed "
            "construction. This closes the additive-IV realization of the A2 "
            "fibre in this X(6,79) route, but does not construct rank 30 and does "
            "not close the separate semistable I3 realization."
        ),
        "theorem": {
            "hypotheses": [
                "c4 and c6 are in the normalized I12+I4+IV logarithmic-derivative model with the three rational fibre locations normalized to infinity, 0, and 1",
                "the discriminant orders are exact and the residual discriminant is nondegenerate",
                "the I4 and IV tangent square classes are both split over Q",
            ],
            "conclusion": (
                "There is no Q-rational surface satisfying all of these hypotheses."
            ),
        },
        "proof_chain": [
            {
                "step": "singular determinant chart",
                "statement": (
                    "On D=0, exact elimination gives a generic resultant with no "
                    "rational roots and a coprime exceptional resultant pair; hence "
                    "there is no Q-rational solution of the surface equations."
                ),
                "source_record_sha256": singular["record_sha256"],
            },
            {
                "step": "open determinant chart",
                "statement": (
                    "On D*kappa!=0, exact QQ Rabinowitsch computations prove R1 "
                    "and R2 lie in the radical of the saturated residual ideal, "
                    "while every residual lies in J=(R1,R2); hence the open locus "
                    "is set-theoretically exactly the known irreducible component."
                ),
                "source_record_sha256": open_component["record_sha256"],
            },
            {
                "step": "split square classes on the unique component",
                "statement": (
                    "The component is parametrized by w over Q. The I4 split square "
                    "class can be positive only for w<-3 or w>3, whereas the IV "
                    "split square class can be positive only for -3<w<0; all "
                    "exceptional parameters are checked separately."
                ),
                "source_record_sha256": split["record_sha256"],
            },
        ],
        "mathematical_consequence": (
            "The additive IV realization cannot supply the split rank-17 K3 seed "
            "of height 79/12, so no section search on that branch is needed."
        ),
        "highest_value_next_action": (
            "Return to the semistable I3 realization of the A2 fibre, but formulate "
            "its surface-and-section incidence using the same sparse Wronskian and "
            "logarithmic-derivative reductions before any new height search."
        ),
        "source_certificates": {
            "singular_determinant": {
                "path": str(singular_path),
                "record_sha256": singular["record_sha256"],
            },
            "open_radical_membership": {
                "path": str(open_path),
                "record_sha256": open_component["record_sha256"],
            },
            "rational_component_split": {
                "path": str(split_path),
                "record_sha256": split["record_sha256"],
            },
        },
        "limitations": [
            "The conclusion concerns the declared normalized additive-IV realization, not every possible strategy for producing rank 30.",
            "The semistable I3 realization of the A2 root fibre remains unresolved.",
            "No elliptic curve over Q with 30 independent points has been produced.",
            "The unconditional global lower bound remains rank E(Q) >= 29.",
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def markdown(payload: dict[str, object]) -> str:
    return "\n".join([
        "# Additive-IV rank-17 branch closure",
        "",
        "```text",
        "SOLVED-A: false",
        "SOLVED-B: false",
        "unconditional global lower bound: rank E(Q) >= 29",
        "```",
        "",
        payload["truth_status"],
        "",
        "## Exact proof chain",
        "",
        "1. **D = 0:** no rational solution of the exact surface equations.",
        "2. **D*kappa != 0:** the saturated open locus is exactly the irreducible rational component `J=(R1,R2)` over the algebraic closure of Q.",
        "3. **On J(Q):** the I4 and IV split square classes require disjoint real intervals, with all exceptional parameters checked separately.",
        "",
        "Therefore:",
        "",
        "> **No nondegenerate Q-rational surface in the normalized**",
        "> **I12 + I4 + IV model satisfies both required split conditions.**",
        "",
        "The additive `IV` realization cannot produce the proposed rank-17 seed.",
        "The separate semistable `I3` realization remains open.",
        "",
        "## Next action",
        "",
        payload["highest_value_next_action"],
        "",
        f"Record SHA-256: `{payload['record_sha256']}`",
        "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--singular", type=Path, required=True)
    parser.add_argument("--open", dest="open_path", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    payload = build(args.singular, args.open_path, args.split)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(markdown(payload), encoding="utf-8")
    if args.compare:
        expected = json.loads(args.compare.read_text(encoding="utf-8"))
        if expected != payload:
            raise AssertionError(f"certificate mismatch: {args.compare}")
    print(json.dumps({
        "classification": payload["classification"],
        "mathematical_consequence": payload["mathematical_consequence"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
