#!/usr/bin/env python3
"""Build a fail-closed checkpoint for the additive-IV rank-17 pipeline."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def load_optional(root: Path, name: str):
    path = root / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def compact_source(data, fields):
    if data is None:
        return {"status": "PENDING"}
    result = {"status": "AVAILABLE"}
    for field in fields:
        result[field] = data.get(field)
    result["source_record_sha256"] = (
        data.get("record_sha256") or data.get("certificate_sha256")
    )
    return result


def build(data_root: Path):
    surfaces11 = load_optional(data_root, "rank17_iv_surfaces_f11_v2.json")
    quotient11 = load_optional(
        data_root, "rank17_iv_surface_geometric_quotient_f11.json"
    )
    endpoint11 = load_optional(
        data_root, "rank17_iv_all_geometric_surfaces_endpoint_f11_v2.json"
    )
    curve11 = load_optional(
        data_root, "rank17_iv_surface_curve_elimination_f11_v2.json"
    )
    surfaces13 = load_optional(data_root, "rank17_iv_surfaces_f13_v2.json")
    quotient13 = load_optional(
        data_root, "rank17_iv_surface_geometric_quotient_f13.json"
    )
    crt = load_optional(
        data_root, "rank17_iv_surface_crt_f11_f13_bounded.json"
    )

    gates = {
        "surface_locus_f11": compact_source(
            surfaces11,
            [
                "prime",
                "visited_parameter_tuples",
                "candidate_count",
                "truth_status",
            ],
        ),
        "surface_geometric_quotient_f11": compact_source(
            quotient11,
            [
                "prime",
                "raw_candidate_count",
                "geometric_surface_count",
                "truth_status",
            ],
        ),
        "section_endpoint_algebras_f11": compact_source(
            endpoint11,
            [
                "classification",
                "geometric_surface_count",
                "expected_sign_orbit_branch_count",
                "observed_branch_count",
                "empty_branch_count",
                "finite_nonempty_branch_count",
                "positive_dimensional_branch_count",
                "incomplete_branch_count",
                "truth_status",
            ],
        ),
        "surface_plane_elimination_f11": compact_source(
            curve11,
            [
                "classification",
                "completed_projection_count",
                "best_plane_equation",
                "truth_status",
            ],
        ),
        "surface_locus_f13": compact_source(
            surfaces13,
            [
                "prime",
                "visited_parameter_tuples",
                "candidate_count",
                "truth_status",
            ],
        ),
        "surface_geometric_quotient_f13": compact_source(
            quotient13,
            [
                "prime",
                "raw_candidate_count",
                "geometric_surface_count",
                "truth_status",
            ],
        ),
        "bounded_crt_f11_f13": compact_source(
            crt,
            [
                "modulus",
                "rational_numerator_denominator_bound",
                "surface_pair_count",
                "exact_q_surface_count",
                "truth_status",
            ],
        ),
    }

    endpoint_classification = (
        endpoint11.get("classification") if endpoint11 else None
    )
    exact_q_surface_count = (
        int(crt.get("exact_q_surface_count", 0)) if crt else 0
    )
    finite_branch_count = (
        int(endpoint11.get("finite_nonempty_branch_count", 0))
        if endpoint11 else 0
    )
    incomplete_branch_count = (
        int(endpoint11.get("incomplete_branch_count", 0))
        if endpoint11 else 0
    )

    if finite_branch_count > 0:
        next_action = (
            "Repeat the endpoint algebra at F13 using the same geometric "
            "coordinates, identify stable quotient degrees and invariant "
            "minimal polynomials, then CRT-reconstruct the section incidence."
        )
    elif exact_q_surface_count > 0:
        next_action = (
            "Reduce each exactly reconstructed Q surface at additional good "
            "primes and solve its saturated section ideal; the surface alone "
            "does not yet realize the height-79/12 generator."
        )
    elif endpoint_classification == "EVERY_GEOMETRIC_F11_SECTION_BRANCH_EMPTY":
        next_action = (
            "Repeat the complete geometric section-branch classification at "
            "F13 and F17 and cover the remaining denominator charts before "
            "treating the F11 emptiness as a serious obstruction."
        )
    elif incomplete_branch_count > 0:
        next_action = (
            "Use the recorded residual degrees and partial bases to split the "
            "incomplete ideals by q=0 versus q!=0 and apply bidirectional square "
            "recursion before another Groebner calculation."
        )
    elif surfaces11 is None:
        next_action = (
            "Complete and publish the independent C/Python enumeration of the "
            "normalized split I12+I4+IV surface locus over F11."
        )
    else:
        next_action = (
            "Complete the geometric quotient and endpoint-reduced section "
            "algebras over F11 without inferring anything from pending jobs."
        )

    payload = {
        "schema_version": 1,
        "checkpoint_id": "rank17_iv_pipeline_checkpoint",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "truth_status": (
            "No rank-30 curve or impossibility theorem is certified.  Every "
            "available gate below is exact only within its declared finite-field, "
            "normalization, chart, and reconstruction bounds."
        ),
        "gates": gates,
        "highest_value_next_action": next_action,
        "global_limitations": [
            "A modular surface is not a characteristic-zero surface.",
            "A modular finite section algebra is not a Q-rational section.",
            "The monic quadratic denominator chart does not by itself cover every projective denominator chart.",
            "The IV branch is only one realization of the A2 fibre; the semistable I3 branch remains separate.",
            "No result in this checkpoint changes the certified rank-29 baseline."
        ],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def markdown(payload):
    lines = [
        "# Rank-17 additive-IV pipeline checkpoint",
        "",
        "```text",
        f"SOLVED-A: {str(payload['solved_a']).lower()}",
        f"SOLVED-B: {str(payload['solved_b']).lower()}",
        (
            "unconditional global lower bound: rank E(Q) >= "
            f"{payload['unconditional_global_rank_lower_bound']}"
        ),
        "```",
        "",
        payload["truth_status"],
        "",
        "## Exact gates",
        "",
    ]
    for name, gate in payload["gates"].items():
        lines.append(f"### `{name}`")
        lines.append("")
        lines.append("```json")
        lines.append(json.dumps(gate, indent=2, sort_keys=True))
        lines.append("```")
        lines.append("")
    lines.extend([
        "## Highest-value next action",
        "",
        payload["highest_value_next_action"],
        "",
        "## Global limitations",
        "",
    ])
    lines.extend(f"- {item}" for item in payload["global_limitations"])
    lines.extend([
        "",
        f"Record SHA-256: `{payload['record_sha256']}`",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=Path("research/data"))
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--markdown-output", type=Path, required=True)
    arguments = parser.parse_args()

    payload = build(arguments.data_root)
    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    arguments.markdown_output.write_text(
        markdown(payload), encoding="utf-8"
    )
    print(json.dumps({
        "available_gate_count": sum(
            gate["status"] == "AVAILABLE"
            for gate in payload["gates"].values()
        ),
        "highest_value_next_action": payload["highest_value_next_action"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
