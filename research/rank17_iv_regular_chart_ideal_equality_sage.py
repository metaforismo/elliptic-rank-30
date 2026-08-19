#!/usr/bin/env sage -python
"""Verify the regular additive-IV chart against its explicit rational component.

The logarithmic-derivative reduction produces a saturated residual ideal I in
F_p[r2,r3,k] after solving the lower linear recurrence away from

    D = -2*k^2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3^2 + 148*r3 - 320.

The proposed regular component is the complete intersection J=(L,F), where L
is linear in r2 and F is the irreducible plane quintic.  This script localizes
both ideals at D with a Rabinowitsch variable h and tests exact two-way ideal
membership using independent Singular standard bases.

A successful finite-field run proves equality only in the declared modular
chart.  It is not by itself a characteristic-zero or rank-30 conclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import time
import traceback
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version

HERE = Path(__file__).resolve().parent
CONSTRUCTION = HERE / "rank17_iv_three_variable_elimination_sage.py"


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_construction_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_three_variable_elimination_sage", CONSTRUCTION
    )
    if spec is None or spec.loader is None:
        raise ImportError(CONSTRUCTION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalize_basis(ring, basis_raw):
    values = [ring(value) for value in basis_raw]
    return [value for value in values if value]


def basis_record(basis) -> dict[str, object]:
    text = "\n".join(str(value) for value in basis) + "\n"
    return {
        "count": len(basis),
        "degrees": [int(value.total_degree()) for value in basis],
        "term_counts": [len(value.dict()) for value in basis],
        "sha256": hashlib.sha256(text.encode()).hexdigest(),
        "text": text,
    }


def compute(prime: int, output: Path) -> dict[str, object]:
    if prime in {2, 3, 5, 11}:
        raise ValueError(("bad prime for this logarithmic chart", prime))

    module = load_construction_module()
    construction = module.construct_residuals(prime)
    source_ring = construction["base"]
    source_r2, source_r3, source_k = source_ring.gens()
    field = GF(prime)

    ring = PolynomialRing(
        field,
        names=("h", "r2", "k", "r3"),
        order="degrevlex",
    )
    h, r2, k, r3 = ring.gens()
    embedding = source_ring.hom([r2, r3, k], ring)

    residuals = [embedding(value) for value in construction["residuals"]]
    D = ring(
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    L = ring(
        9*k**3 + 117*k**2*r3 - 194*k**2
        + 386*k*r3 + 436*k + 5400*r2
        - 950*r3**2 + 4072*r3 + 5824
    )
    F = ring(
        81*k**5 + 1053*k**4*r3 - 1098*k**4
        + 13248*k**3*r3 - 11048*k**3
        + 9000*k**2*r3**2 + 116688*k**2*r3 + 115872*k**2
        + 12000*k*r3**2 + 87168*k*r3 + 112512*k
        + 40000*r3**3 + 384000*r3**2 + 979968*r3 + 813056
    )

    inverse_equation = h*D - 1
    residual_generators = residuals + [inverse_equation]
    component_generators = [L, F, inverse_equation]

    initial = {
        "schema_version": 1,
        "status": "standard_bases_started",
        "truth_status": (
            f"exact ideal comparison on the D!=0 logarithmic chart over F_{prime}; "
            "no characteristic-zero, section, or rank-30 conclusion"
        ),
        "prime": prime,
        "sage_version": str(sage_version),
        "ring_variables": [str(value) for value in ring.gens()],
        "residual_equation_count": len(residuals),
        "residual_equation_degrees": [
            int(value.total_degree()) for value in residuals
        ],
        "candidate_component": {
            "linear_equation": str(L),
            "quintic_equation": str(F),
            "localizing_polynomial": str(D),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")

    slimgb = singular_function("slimgb")
    reduce_function = singular_function("reduce")
    dim_function = singular_function("dim")

    started = time.time()
    residual_basis_raw = slimgb(ring.ideal(residual_generators))
    residual_elapsed = time.time() - started
    residual_basis = normalize_basis(ring, residual_basis_raw)

    started = time.time()
    component_basis_raw = slimgb(ring.ideal(component_generators))
    component_elapsed = time.time() - started
    component_basis = normalize_basis(ring, component_basis_raw)

    residual_record = basis_record(residual_basis)
    component_record = basis_record(component_basis)
    (output.parent / "regular-residual-basis.txt").write_text(
        residual_record.pop("text")
    )
    (output.parent / "regular-component-basis.txt").write_text(
        component_record.pop("text")
    )

    component_in_residual = []
    for generator in component_generators:
        remainder = ring(reduce_function(generator, residual_basis_raw))
        component_in_residual.append({
            "generator": str(generator),
            "remainder": str(remainder),
            "is_zero": not bool(remainder),
        })

    residual_in_component = []
    for index, generator in enumerate(residual_generators):
        remainder = ring(reduce_function(generator, component_basis_raw))
        residual_in_component.append({
            "generator_index": index,
            "remainder": str(remainder),
            "is_zero": not bool(remainder),
        })

    exact_equality = (
        all(item["is_zero"] for item in component_in_residual)
        and all(item["is_zero"] for item in residual_in_component)
    )

    result = dict(initial)
    result.update({
        "status": "completed",
        "residual_basis_elapsed_seconds": residual_elapsed,
        "component_basis_elapsed_seconds": component_elapsed,
        "residual_basis": residual_record,
        "component_basis": component_record,
        "residual_krull_dimension": int(dim_function(residual_basis_raw)),
        "component_krull_dimension": int(dim_function(component_basis_raw)),
        "candidate_generators_reduce_to_zero_mod_residual": component_in_residual,
        "residual_generators_reduce_to_zero_mod_candidate": residual_in_component,
        "localized_ideals_exactly_equal": exact_equality,
        "limitations": [
            "The determinant-zero recurrence chart is excluded and must be analysed separately.",
            "Finite-field equality at one prime is not a characteristic-zero equality certificate.",
            "No split-square condition or height-79/12 section is imposed here.",
            "The unconditional global rank lower bound remains 29."
        ],
    })
    result["record_sha256"] = canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        result = compute(arguments.prime, arguments.output)
        print(json.dumps({
            "status": result["status"],
            "prime": result["prime"],
            "localized_ideals_exactly_equal": result[
                "localized_ideals_exactly_equal"
            ],
            "residual_krull_dimension": result["residual_krull_dimension"],
            "record_sha256": result["record_sha256"],
        }, indent=2, sort_keys=True), flush=True)
    except Exception as exc:
        error = {
            "schema_version": 1,
            "status": "error",
            "truth_status": "implementation or computation failure; no mathematical conclusion",
            "prime": arguments.prime,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
        }
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(error, indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(error, indent=2, sort_keys=True), flush=True)
        raise


if __name__ == "__main__":
    main()
