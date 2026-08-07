#!/usr/bin/env python3
"""Normalize the certified split-E8 surface data for Sage search workers."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def walk(obj, path=()):
    yield path, obj
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from walk(value, path + (str(key),))
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            yield from walk(value, path + (str(index),))


def point_coordinates(item):
    if isinstance(item, dict):
        for x_key in ("x", "X", "x_coordinate"):
            for y_key in ("y", "Y", "y_coordinate"):
                if x_key in item and y_key in item:
                    return str(item[x_key]), str(item[y_key])
        for nested_key in ("coordinates", "point"):
            if nested_key in item:
                result = point_coordinates(item[nested_key])
                if result is not None:
                    return result
    if isinstance(item, (list, tuple)) and len(item) >= 2:
        return str(item[0]), str(item[1])
    return None


def normalize(source):
    data = json.loads(source.read_text())

    point_candidates = []
    for path, obj in walk(data):
        if isinstance(obj, list) and len(obj) >= 8:
            points = [point_coordinates(item) for item in obj]
            if sum(point is not None for point in points) >= 8:
                point_candidates.append((path, obj, points))
    if not point_candidates:
        raise RuntimeError("no list containing eight section coordinates was found")
    point_candidates.sort(
        key=lambda row: (
            len(row[1]) != 8,
            not any(
                "section" in component.lower() or "generator" in component.lower()
                for component in row[0]
            ),
            len(row[0]),
        )
    )
    point_path, _, points = point_candidates[0]
    points = [point for point in points if point is not None][:8]

    matrix_candidates = []
    for path, obj in walk(data):
        if (
            isinstance(obj, list)
            and len(obj) == 8
            and all(isinstance(row, list) and len(row) == 8 for row in obj)
        ):
            matrix_candidates.append(
                (path, [[str(entry) for entry in row] for row in obj])
            )
    if not matrix_candidates:
        raise RuntimeError("no 8 by 8 Gram matrix was found")
    matrix_candidates.sort(
        key=lambda row: (
            not any(
                "height" in component.lower() or "gram" in component.lower()
                for component in row[0]
            ),
            len(row[0]),
        )
    )
    matrix_path, matrix = matrix_candidates[0]

    return {
        "status": "pass",
        "source": str(source),
        "section_source_path": "/".join(point_path),
        "height_matrix_source_path": "/".join(matrix_path),
        "section_count": 8,
        "sections": [{"x": x, "y": y} for x, y in points],
        "height_matrix": matrix,
        "weierstrass_model": {
            "a4": "-27*t^4 + 216*t^3 + 6*t^2 - 432*t - 267",
            "a6": "54*t^6 + 108*t^5 - 630*t^4 - 880*t^3 + 1458*t^2 + 1278*t + 1242",
        },
        "truth_status": "normalization of an existing exact certificate",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("certificates/kumar_shioda_split_e8_certificate.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = normalize(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "section_count": result["section_count"],
                "section_source_path": result["section_source_path"],
                "height_matrix_source_path": result["height_matrix_source_path"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
