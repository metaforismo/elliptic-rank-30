#!/usr/bin/env python3
"""Build or verify the canonical SHA-256 manifest of tracked repository files.

The manifest deliberately excludes itself and Git internals.  Every other
tracked file, including generated proof artifacts, the Lean dependency graph,
and the canonical paper PDF, is covered once present on the branch.
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELF = "MANIFEST.sha256"


def tracked_paths() -> list[str]:
    output = subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=ROOT
    )
    paths = [item.decode("utf-8") for item in output.split(b"\0") if item]
    return sorted(path for path in paths if path != SELF)


def render_manifest() -> str:
    rows: list[str] = []
    for relative in tracked_paths():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"tracked path is not a regular file: {relative}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{digest}  {relative}")
    return "\n".join(rows) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="write the generated manifest to this path",
    )
    parser.add_argument(
        "--check",
        type=Path,
        help="compare the generated manifest with this existing file",
    )
    args = parser.parse_args()

    if args.output and args.check:
        parser.error("--output and --check are mutually exclusive")

    generated = render_manifest()
    if args.check:
        expected = (ROOT / args.check).read_text()
        if generated != expected:
            generated_lines = generated.splitlines()
            expected_lines = expected.splitlines()
            limit = max(len(generated_lines), len(expected_lines))
            for index in range(limit):
                actual = generated_lines[index] if index < len(generated_lines) else "<missing>"
                wanted = expected_lines[index] if index < len(expected_lines) else "<missing>"
                if actual != wanted:
                    raise SystemExit(
                        f"manifest mismatch at line {index + 1}:\n"
                        f"expected: {wanted}\nactual:   {actual}"
                    )
            raise SystemExit("manifest mismatch")
        print(f"verified {len(generated.splitlines())} tracked files")
        return 0

    if args.output:
        destination = ROOT / args.output
        destination.write_text(generated)
        print(f"wrote {destination} with {len(generated.splitlines())} entries")
    else:
        print(generated, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
