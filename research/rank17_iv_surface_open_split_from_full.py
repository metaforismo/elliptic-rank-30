#!/usr/bin/env python3
"""Extract the exact split open IV locus from a full finite-field census.

Input records use the e0=1 jet-sign normalization and include every exact
I12+I4+IV surface before tangent-square filtering.  This script retains only
surfaces for which both tangent targets split over F_p and for which the
residual discriminant quartic is degree four, squarefree, and coprime to c4.
Thus the accepted records have exactly

    I12 + I4 + IV + 4 I1

in the declared normalized chart.  This remains finite-field evidence only.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESIDUAL_SCRIPT = HERE / "rank17_iv_surface_residual_geometry.py"


def load_residual_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_surface_residual_geometry", RESIDUAL_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(RESIDUAL_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build(full_path: Path) -> dict[str, object]:
    residual_module = load_residual_module()
    full = json.loads(full_path.read_text(encoding="utf-8"))
    prime = int(full["prime"])
    split_records = [
        record for record in full["records"]
        if record["split_tangent_checks"]["i4_split"]
        and record["split_tangent_checks"]["iv_split"]
    ]
    analyses = [
        residual_module.analyze_candidate(record, prime)
        for record in split_records
    ]
    accepted_hashes = {
        analysis["representative_record_sha256"]
        for analysis in analyses
        if analysis["accepted_exact_fibre_configuration"]
    }
    representatives = [
        record for record in split_records
        if record["record_sha256"] in accepted_hashes
    ]
    degenerate = [
        analysis for analysis in analyses
        if not analysis["accepted_exact_fibre_configuration"]
    ]
    payload = {
        "schema_version": 1,
        "certificate_id": f"rank17_iv_open_split_surface_locus_f{prime}",
        "truth_status": (
            f"EXHAUSTIVE split open I12+I4+IV+4I1 surface locus over F_{prime} "
            "inside the e0=1 normalized chart; no characteristic-zero, section, "
            "or rank-30 conclusion"
        ),
        "prime": prime,
        "normalization": full["normalization"],
        "source_full_locus_sha256": full["record_sha256"],
        "full_surface_count": full["candidate_count"],
        "both_tangents_split_count": len(split_records),
        "accepted_exact_configuration_count": len(representatives),
        "rejected_residual_degeneration_count": len(degenerate),
        "exact_configuration": "I12 + I4 + IV + 4 I1",
        "analyses": analyses,
        "representatives": representatives,
        "limitations": [
            "The census is confined to the declared normalized affine chart.",
            "A finite-field point need not lift to characteristic zero.",
            "No height-79/12 section is imposed.",
            "The unconditional global rank lower bound remains 29."
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build(arguments.full)
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
        "prime": payload["prime"],
        "full_surface_count": payload["full_surface_count"],
        "both_tangents_split_count": payload["both_tangents_split_count"],
        "accepted_exact_configuration_count": payload["accepted_exact_configuration_count"],
        "rejected_residual_degeneration_count": payload["rejected_residual_degeneration_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
