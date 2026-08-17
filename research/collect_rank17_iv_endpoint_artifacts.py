#!/usr/bin/env python3
"""Fail-closed collector for endpoint-reduced IV section artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_artifact_coordinates(name: str) -> dict[str, int | None]:
    parts = name.rsplit("-", 3)
    result: dict[str, int | None] = {
        "candidate_index": None,
        "i4_root": None,
        "iv_root": None,
    }
    if len(parts) != 4:
        return result
    try:
        result.update({
            "candidate_index": int(parts[1]),
            "i4_root": int(parts[2]),
            "iv_root": int(parts[3]),
        })
    except ValueError:
        pass
    return result


def collect_branch_records(root: Path) -> list[dict[str, object]]:
    records = []
    if not root.exists():
        return records
    for result_path in sorted(root.rglob("result.json")):
        directory = result_path.parent
        result = json.loads(result_path.read_text(encoding="utf-8"))
        coordinates = parse_artifact_coordinates(directory.name)
        for key in coordinates:
            if result.get(key) is not None:
                coordinates[key] = int(result[key])
        basis_path = directory / "groebner-basis.txt"
        basis_text = basis_path.read_text() if basis_path.exists() else None
        exit_path = directory / "exit-status.txt"
        exit_status = None
        if exit_path.exists() and exit_path.read_text().strip():
            exit_status = int(exit_path.read_text().strip())
        records.append({
            **coordinates,
            "artifact_directory": directory.name,
            "status": result.get("status", "missing_status"),
            "truth_status": result.get("truth_status"),
            "exit_status": exit_status,
            "krull_dimension": result.get("krull_dimension"),
            "quotient_dimension": result.get("quotient_dimension"),
            "primitive_etale_coordinate": result.get(
                "primitive_etale_coordinate"
            ),
            "multiplication_coordinates": result.get(
                "multiplication_coordinates", []
            ),
            "groebner_basis_is_unit": (
                basis_text.strip() == "1" if basis_text is not None else None
            ),
            "result_sha256": file_sha256(result_path),
            "basis_sha256": file_sha256(basis_path) if basis_path.exists() else None,
            "error": result.get("error"),
        })
    records.sort(key=lambda item: (
        item["candidate_index"] if item["candidate_index"] is not None else -1,
        item["i4_root"] if item["i4_root"] is not None else -1,
        item["iv_root"] if item["iv_root"] is not None else -1,
    ))
    return records


def build(
    *,
    quotient_path: Path,
    artifacts_root: Path,
    prepare_result: str,
    solve_result: str,
) -> dict[str, object]:
    quotient = json.loads(quotient_path.read_text(encoding="utf-8"))
    records = collect_branch_records(artifacts_root)
    surface_count = int(quotient["geometric_surface_count"])
    expected = 2 * surface_count

    empty = [
        item for item in records
        if item["groebner_basis_is_unit"] is True
        or (
            isinstance(item["krull_dimension"], int)
            and item["krull_dimension"] < 0
        )
        or item["quotient_dimension"] == 0
    ]
    finite = [
        item for item in records
        if item["krull_dimension"] == 0
        and isinstance(item["quotient_dimension"], int)
        and item["quotient_dimension"] > 0
    ]
    positive = [
        item for item in records
        if isinstance(item["krull_dimension"], int)
        and item["krull_dimension"] > 0
    ]
    incomplete = [
        item for item in records
        if item["status"] != "completed"
        or item["exit_status"] not in (0, None)
    ]

    if surface_count == 0:
        classification = "NO_GEOMETRIC_SPLIT_IV_SURFACE_OVER_F11"
    elif finite:
        classification = "AT_LEAST_ONE_FINITE_SECTION_ALGEBRA_FOUND"
    elif positive:
        classification = "AT_LEAST_ONE_POSITIVE_DIMENSIONAL_SECTION_LOCUS"
    elif len(records) == expected and len(empty) == expected and not incomplete:
        classification = "EVERY_GEOMETRIC_F11_SECTION_BRANCH_EMPTY"
    elif surface_count > 128:
        classification = "TOO_MANY_SURFACES_FOR_SINGLE_MATRIX_RUN"
    else:
        classification = "ALL_SURFACE_BRANCH_CLASSIFICATION_INCOMPLETE"

    payload: dict[str, object] = {
        "schema_version": 2,
        "certificate_id": "rank17_iv_all_geometric_surfaces_endpoint_f11_v2",
        "truth_status": (
            "EXACT F11 endpoint-reduced section classification where completed; "
            "timeouts and failures remain unresolved; no p-adic, characteristic-zero, "
            "rational-section, or rank-30 conclusion"
        ),
        "geometric_quotient_sha256": quotient["record_sha256"],
        "geometric_surface_count": surface_count,
        "expected_sign_orbit_branch_count": expected,
        "observed_branch_count": len(records),
        "classification": classification,
        "empty_branch_count": len(empty),
        "finite_nonempty_branch_count": len(finite),
        "positive_dimensional_branch_count": len(positive),
        "incomplete_branch_count": len(incomplete),
        "workflow_needs": {
            "prepare": prepare_result,
            "solve": solve_result,
        },
        "branches": records,
        "limitations": [
            "This covers the normalized monic-D chart over F11 only.",
            "Section negation is quotiented by fixing one I4 root and testing both relative IV signs.",
            "A finite modular algebra is not a Q-rational point.",
            "A timeout or failed branch remains mathematically open.",
            "The unconditional global rank lower bound remains 29.",
        ],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quotient", type=Path, required=True)
    parser.add_argument("--artifacts", type=Path, required=True)
    parser.add_argument("--prepare-result", default="unknown")
    parser.add_argument("--solve-result", default="unknown")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    payload = build(
        quotient_path=arguments.quotient,
        artifacts_root=arguments.artifacts,
        prepare_result=arguments.prepare_result,
        solve_result=arguments.solve_result,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "classification": payload["classification"],
        "surfaces": payload["geometric_surface_count"],
        "observed": payload["observed_branch_count"],
        "empty": payload["empty_branch_count"],
        "finite": payload["finite_nonempty_branch_count"],
        "positive": payload["positive_dimensional_branch_count"],
        "incomplete": payload["incomplete_branch_count"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
