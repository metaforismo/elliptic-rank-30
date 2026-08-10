#!/usr/bin/env sage-python
"""Replay and serialize the deterministic neighbor chain to the X(6,79) lattice.

This utility deliberately derives the discovery program from the exact workflow
that produced the candidate, patches only reproducibility/output concerns, and
requires the previously frozen target Gram-matrix hash.  It does not assume the
mathematical conclusion merely because the discovery run succeeded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import textwrap
from pathlib import Path

EXPECTED_TARGET_HASH = "620a5e06473684d3e8015c0172f63c09c901e742ec02e77ba0aa35a923aa0295"
DEFAULT_SOURCE = Path(".github/workflows/probe-rank29-height-lattice.yml")


def extract_discovery_script(source: Path) -> str:
    workflow = source.read_text()
    start_marker = "cat > /tmp/rank17_root_elimination.py <<'PY'\n"
    end_marker = "\n          PY\n          sage -python /tmp/rank17_root_elimination.py"
    if start_marker not in workflow or end_marker not in workflow:
        raise RuntimeError("could not locate the embedded root-elimination program")
    return textwrap.dedent(workflow.split(start_marker, 1)[1].split(end_marker, 1)[0])


def patch_script(script: str, output: Path) -> str:
    old_import = (
        "from sage.all import Genus, IntegralLattice, Matrix, QuadraticForm, ZZ, "
        "block_diagonal_matrix, matrix, pari, random_seed, version"
    )
    new_import = (
        "from sage.all import Genus, IntegralLattice, Matrix, QuadraticForm, ZZ, "
        "block_diagonal_matrix, matrix, pari, version\n"
        "from sage.misc.randstate import set_random_seed"
    )
    if old_import not in script:
        raise RuntimeError("the discovery program no longer has the expected Sage import")
    script = script.replace(old_import, new_import, 1)
    script = script.replace("random_seed(RANDOM_SEED)", "set_random_seed(RANDOM_SEED)", 1)

    old_output = "OUT = Path('/tmp/rank17-root-elimination')"
    new_output = f"OUT = Path({str(output)!r})"
    if old_output not in script:
        raise RuntimeError("the discovery output assignment changed unexpectedly")
    script = script.replace(old_output, new_output, 1)

    result_needle = """              result.update({
                  'status': 'target_found' if verified_targets else 'bounded_search_completed',
"""
    result_replacement = """              state_by_hash = {state['hash']: state for state in archive}
              target_neighbor_chains = []
              for target_state in targets:
                  chain = []
                  cursor = target_state
                  while cursor is not None:
                      chain.append(public_record(cursor, index_by_hash))
                      cursor = state_by_hash.get(cursor['parent']) if cursor['parent'] is not None else None
                  chain.reverse()
                  target_neighbor_chains.append(chain)

              result.update({
                  'status': 'target_found' if verified_targets else 'bounded_search_completed',
"""
    if result_needle not in script:
        raise RuntimeError("could not locate the final result construction")
    script = script.replace(result_needle, result_replacement, 1)

    field_needle = """                  'verified_targets': verified_targets,
                  'target_found': bool(verified_targets),
"""
    field_replacement = """                  'verified_targets': verified_targets,
                  'target_neighbor_chains': target_neighbor_chains,
                  'target_found': bool(verified_targets),
"""
    if field_needle not in script:
        raise RuntimeError("could not add neighbor chains to the discovery record")
    script = script.replace(field_needle, field_replacement, 1)
    return script


def canonical_matrix_hash(matrix_data: list[list[int]]) -> str:
    raw = json.dumps(matrix_data, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    script = patch_script(extract_discovery_script(args.source), output)
    replay_script = output / "replayed-discovery.py"
    replay_script.write_text(script)

    completed = subprocess.run(
        [sys.executable, str(replay_script)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (output / "replay.log").write_text(completed.stdout)
    print(completed.stdout, end="")
    if completed.returncode:
        raise SystemExit(f"discovery replay failed with exit code {completed.returncode}")

    record_path = output / "rank17-root-elimination.json"
    record = json.loads(record_path.read_text())
    if record.get("status") != "target_found":
        raise SystemExit(f"replay did not reach the frozen target: {record.get('status')}")
    verified = record.get("verified_targets") or []
    hashes = [entry.get("hash") for entry in verified]
    if EXPECTED_TARGET_HASH not in hashes:
        raise SystemExit(f"replay reached different target hashes: {hashes}")

    chains = record.get("target_neighbor_chains") or []
    matching = [chain for chain in chains if chain and chain[-1].get("hash") == EXPECTED_TARGET_HASH]
    if len(matching) != 1:
        raise SystemExit(f"expected one serialized chain to the frozen target, found {len(matching)}")
    chain = matching[0]
    if any(chain[index].get("parent_hash") != chain[index - 1].get("hash") for index in range(1, len(chain))):
        raise SystemExit("serialized parent links do not form a chain")
    if canonical_matrix_hash(chain[-1]["gram"]) != EXPECTED_TARGET_HASH:
        raise SystemExit("terminal Gram matrix does not have the frozen canonical hash")

    exact_chain_path = output / "exact-neighbor-chain.json"
    exact_chain_path.write_text(json.dumps(chain, indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "completed",
        "source_workflow": str(args.source),
        "source_workflow_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "replay_script_sha256": hashlib.sha256(replay_script.read_bytes()).hexdigest(),
        "chain_length_including_seed": len(chain),
        "neighbor_step_count": len(chain) - 1,
        "root_counts": [entry["root_count"] for entry in chain],
        "norm4_counts": [entry["norm4_count"] for entry in chain],
        "neighbor_primes": [entry["move"]["prime"] for entry in chain[1:]],
        "target_hash": EXPECTED_TARGET_HASH,
        "target_gram_sha256": canonical_matrix_hash(chain[-1]["gram"]),
        "chain_sha256": hashlib.sha256(exact_chain_path.read_bytes()).hexdigest(),
    }
    (output / "chain-summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
