#!/usr/bin/env sage -python
"""Exact radical-membership test for the open additive-IV surface locus.

The determinant-nonzero logarithmic-derivative chart has a known rational
component J=(R1,R2). The residual surface equations are already known to lie
in J. To exclude every additional component and isolated point on the open
chart, this script tests over QQ that R1 and R2 belong to the radical of the
saturated residual ideal I by the Rabinowitsch criterion:

    1 in I + (u*R1-1),   1 in I + (u*R2-1).

If both exact unit-ideal tests pass, V(I)=V(J) set-theoretically on the declared
open chart. No section or rank-30 claim follows without the separate split and
section arguments.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "rank17_iv_three_variable_open_elimination_sage.py"


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_constructor():
    source = IMPLEMENTATION.read_text(encoding="utf-8")
    replacements = [
        (
            "from sage.all import GF, PolynomialRing",
            "from sage.all import GF, PolynomialRing, QQ",
        ),
        (
            "    if prime in BAD_PRIMES:\n        raise ValueError((\"bad prime for this coordinate chart\", prime))",
            "    if prime != 0 and prime in BAD_PRIMES:\n        raise ValueError((\"bad prime for this coordinate chart\", prime))",
        ),
        (
            "    field = GF(prime)",
            "    field = QQ if prime == 0 else GF(prime)",
        ),
        (
            """        while True:\n            quotient, remainder = numerator.quo_rem(expected_determinant)\n            if remainder:\n                break\n            numerator = quotient\n            determinant_power += 1\n""",
            """        while True:\n            candidate = fraction(numerator) / fraction(expected_determinant)\n            if candidate.denominator() != 1:\n                break\n            numerator = base(candidate.numerator())\n            determinant_power += 1\n""",
        ),
    ]
    for old, new in replacements:
        if source.count(old) != 1:
            raise AssertionError(("constructor compatibility site changed", old[:80]))
        source = source.replace(old, new)
    namespace = {
        "__name__": "rank17_iv_open_constructor_QQ",
        "__file__": str(IMPLEMENTATION),
    }
    exec(compile(source, str(IMPLEMENTATION), "exec"), namespace)
    return namespace["construct_residuals"]


def polynomial_record(polynomial) -> dict[str, object]:
    return {
        "expression": str(polynomial),
        "total_degree": int(polynomial.total_degree()),
        "term_count": len(polynomial.dict()),
        "sha256": hashlib.sha256(str(polynomial).encode()).hexdigest(),
    }


def compute(output: Path) -> dict[str, object]:
    construct_residuals = load_constructor()
    construction = construct_residuals(0)
    base = construction["base"]
    r2b, r3b, kb = base.gens()

    R1b = base(
        -8496*kb**2*r3b + 4374*kb**2*r2b - 938*kb**2
        - 2025*kb*r3b**2 + 20552*kb*r3b + 34992*kb*r2b
        + 23236*kb - 16425*r3b**3 + 29950*r3b**2
        + 72900*r3b*r2b - 50456*r3b - 54216*r2b - 131648
    )
    R2b = base(
        2832*kb**3 + 18954*kb**2*r2b - 65110*kb**2
        - 8775*kb*r3b**2 + 210520*kb*r3b + 151632*kb*r2b
        + 237884*kb - 71175*r3b**3 - 169150*r3b**2
        + 315900*r3b*r2b + 1062680*r3b + 1464264*r2b
        + 1262144
    )

    component_ideal = base.ideal([R1b, R2b])
    residual_remainders = [
        component_ideal.reduce(polynomial)
        for polynomial in construction["residuals"]
    ]
    if any(residual_remainders):
        raise AssertionError("a residual equation is not in the known component ideal")

    resultant = R1b.resultant(R2b, r2b)
    plane = PolynomialRing(QQ, names=("k", "r3"), order="lex")
    kk, yy = plane.gens()
    resultant_plane = plane(str(resultant)).monic()
    factors = list(resultant_plane.factor())
    if len(factors) != 1 or int(factors[0][1]) != 1:
        raise AssertionError("the component resultant is not irreducible over QQ")

    ring = PolynomialRing(
        QQ,
        names=("u", "h", "r2", "r3", "k"),
        order="degrevlex",
    )
    u, h, r2, r3, k = ring.gens()
    embedding = base.hom([r2, r3, k], ring)
    open_product = embedding(construction["determinant"])*embedding(
        construction["kappa_numerator"]
    )
    residuals = [embedding(value) for value in construction["residuals"]]
    saturation = h*open_product - 1
    targets = {
        "R1": embedding(R1b),
        "R2": embedding(R2b),
    }

    slimgb = singular_function("slimgb")
    target_results = {}
    output.parent.mkdir(parents=True, exist_ok=True)
    for name, target in targets.items():
        initial = ring.ideal(residuals + [saturation, u*target - 1])
        started = time.time()
        basis_raw = slimgb(initial)
        elapsed = time.time() - started
        basis = [ring(value) for value in basis_raw]
        unit = any(value == 1 for value in basis)
        basis_text = "\n".join(str(value) for value in basis) + "\n"
        (output.parent / f"radical-membership-{name}.txt").write_text(
            basis_text, encoding="utf-8"
        )
        if not unit:
            raise AssertionError(("radical membership failed", name))
        target_results[name] = {
            "rabinowitsch_unit_ideal": True,
            "elapsed_seconds": elapsed,
            "basis_count": len(basis),
            "basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        }

    result: dict[str, object] = {
        "schema_version": 1,
        "certificate_id": "rank17_iv_open_component_radical_membership",
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "truth_status": (
            "CERTIFIED over QQ: on the declared determinant- and discriminant-"
            "nonzero logarithmic-derivative chart, every solution of the exact "
            "additive-IV surface equations lies on the known irreducible rational "
            "component J=(R1,R2). This is not by itself a rank-30 construction."
        ),
        "sage_version": str(sage_version),
        "residual_count": len(residuals),
        "residuals_contained_in_J": True,
        "open_saturation_product": str(open_product),
        "component": {
            "R1": polynomial_record(R1b),
            "R2": polynomial_record(R2b),
            "plane_resultant": polynomial_record(resultant_plane),
            "plane_resultant_irreducible_over_QQ": True,
        },
        "radical_membership": target_results,
        "set_theoretic_consequence": (
            "V(I_open)=V(R1,R2)_open over the algebraic closure of Q."
        ),
        "limitations": [
            "The statement is set-theoretic; it does not claim equality of nonreduced schemes.",
            "The determinant-zero chart is handled by a separate exact rational obstruction.",
            "The known component still must satisfy the two split-fibre square conditions; a separate certificate proves they are incompatible over Q.",
            "No height-79/12 section or rank-30 elliptic curve is constructed.",
            "The semistable I3 realization remains separate.",
        ],
    }
    result["record_sha256"] = canonical_hash(result)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = compute(arguments.output)
    print(json.dumps({
        "radical_membership": result["radical_membership"],
        "set_theoretic_consequence": result["set_theoretic_consequence"],
        "record_sha256": result["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
