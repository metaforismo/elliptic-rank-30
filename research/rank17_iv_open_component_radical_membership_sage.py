#!/usr/bin/env sage -python
"""Exact radical-membership certificate for the open additive-IV surface locus.

This file deliberately constructs the logarithmic-derivative residual ideal
*directly over QQ*.  It does not patch or execute a finite-field wrapper.

Write

    R = r0 + r1*t + r2*t^2 + r3*t^3 + t^4,
    P = t*R,
    Q = P' + 3*R,
    K = k - 12*t,
    M = m0 + m1*t + m2*t^2 + 80*t^3,
    N = Q^2 + M*P,

and

    T = N*K*((t-1)*Q - 2*P)
        - 3*P*(t-1)*(N'*K - 2*N*K').

The differential equation is

    Q*T*K - 2*P*T'*K + 8*P*T*K' - (t-1)*N^2*K^2 = 0.

The coefficients of its quotient by P solve m2,m1 and then m0,r1,r0.
After removing powers of the lower determinant D, the remaining exact
polynomials in QQ[r2,r3,k] generate the regular-chart residual ideal.
The coefficient of t^16 in

    (t-1)^2*N^3*K^2 - T^2

gives, up to (-12)^8, the discriminant scalar kappa.  We localize at D*kappa.

For the known rational component J=(R1,R2), the script proves:

  1. every residual generator lies in J;
  2. R1 and R2 lie in the radical of the localized residual ideal I,
     using the Rabinowitsch tests

         1 in I + (u*R1 - 1),
         1 in I + (u*R2 - 1).

Consequently V(I)=V(J) on the declared determinant- and discriminant-open
chart over the algebraic closure of Q.  This is a surface-locus statement;
it constructs no Mordell-Weil section and no rank-30 curve.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ, vector
from sage.libs.singular.function_factory import singular_function
from sage.version import version as sage_version


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def exact_remove_factor(polynomial, factor):
    """Remove the largest exact nonnegative power of factor."""
    value = polynomial
    power = 0
    while value:
        quotient, remainder = value.quo_rem(factor)
        if remainder:
            break
        value = quotient
        power += 1
    return value, power


def construct_open_model_over_qq() -> dict[str, object]:
    base = PolynomialRing(
        QQ,
        names=("r2", "r3", "k"),
        order="degrevlex",
    )
    r2, r3, k = base.gens()
    fraction = base.fraction_field()
    time_ring = PolynomialRing(fraction, "t")
    t = time_ring.gen()
    s = t - 1

    m2 = fraction(8 * (-k + 2*r3 + 4))
    m1 = fraction(-(
        k**2 + 14*k*r3 + 4*k + 24*r2
        + 13*r3**2 - 56*r3 - 80
    )) / 3

    def qbar(*, r0_value=0, r1_value=0, m0_value=0):
        r0 = fraction(r0_value)
        r1 = fraction(r1_value)
        m0 = fraction(m0_value)
        R = r0 + r1*t + r2*t**2 + r3*t**3 + t**4
        P = t*R
        Q = P.derivative() + 3*R
        K = k - 12*t
        M = m0 + m1*t + m2*t**2 + 80*t**3
        N = Q**2 + M*P
        T = (
            N*K*(s*Q - 2*P)
            - 3*P*s*(N.derivative()*K - 2*N*K.derivative())
        )
        E = (
            Q*T*K - 2*P*T.derivative()*K
            + 8*P*T*K.derivative() - s*N**2*K**2
        )
        quotient, remainder = E.quo_rem(P)
        if remainder:
            raise AssertionError("the differential equation lost its forced P factor")
        if quotient.degree() != 13:
            raise AssertionError(("unexpected quotient degree", quotient.degree()))
        return quotient

    degrees = (11, 10, 9)
    zero = qbar()
    base_values = vector(fraction, [zero[degree] for degree in degrees])
    columns = []
    for keyword in ("m0_value", "r1_value", "r0_value"):
        trial = qbar(**{keyword: 1})
        columns.append(vector(
            fraction,
            [trial[degree] - zero[degree] for degree in degrees],
        ))
    matrix = Matrix(
        fraction,
        3,
        3,
        lambda row, column: columns[column][row],
    )
    determinant = matrix.det()
    if not determinant:
        raise AssertionError("the lower recurrence determinant vanished identically")

    D = base(
        -2*k**2 - 37*k*r3 + 52*k + 108*r2
        - 116*r3**2 + 148*r3 - 320
    )
    ratio = fraction(determinant) / fraction(D)
    if ratio.numerator().total_degree() != 0 or ratio.denominator().total_degree() != 0:
        raise AssertionError(("unexpected determinant factor", determinant, ratio))

    solution = matrix.solve_right(-base_values)
    m0, r1, r0 = solution
    final = qbar(r0_value=r0, r1_value=r1, m0_value=m0)
    for degree in (13, 12, 11, 10, 9):
        if final[degree]:
            raise AssertionError(("a solved coefficient survived", degree, final[degree]))

    residuals = []
    residual_metadata = []
    seen = set()
    for degree in range(8, -1, -1):
        value = final[degree]
        numerator = base(value.numerator())
        denominator = base(value.denominator())
        if not numerator:
            continue
        numerator, removed_power = exact_remove_factor(numerator, D)
        numerator = numerator.monic()
        key = str(numerator)
        if key in seen:
            continue
        seen.add(key)
        residuals.append(numerator)
        residual_metadata.append({
            "differential_degree": degree,
            "total_degree": int(numerator.total_degree()),
            "degree_r2": int(numerator.degree(r2)),
            "degree_r3": int(numerator.degree(r3)),
            "degree_k": int(numerator.degree(k)),
            "term_count": len(numerator.dict()),
            "removed_determinant_power": removed_power,
            "original_denominator": str(denominator),
            "sha256": hashlib.sha256(str(numerator).encode()).hexdigest(),
        })
    if not residuals:
        raise AssertionError("no residual polynomial survived")

    R = r0 + r1*t + r2*t**2 + r3*t**3 + t**4
    P = t*R
    Q = P.derivative() + 3*R
    K = k - 12*t
    M = m0 + m1*t + m2*t**2 + 80*t**3
    N = Q**2 + M*P
    T = (
        N*K*(s*Q - 2*P)
        - 3*P*s*(N.derivative()*K - 2*N*K.derivative())
    )
    J = s**2*N**3*K**2 - T**2
    kappa_coefficient = J[16] / ((-12)**8)
    kappa_numerator = base(kappa_coefficient.numerator())
    kappa_denominator = base(kappa_coefficient.denominator())
    if not kappa_numerator:
        raise AssertionError("the open discriminant scalar vanished identically")
    kappa_numerator, kappa_removed_power = exact_remove_factor(
        kappa_numerator, D
    )
    kappa_numerator = kappa_numerator.monic()

    # Replay the defining identity at the generic level.  The lower solution
    # and residuals ensure the differential equation; this coefficient is the
    # exact leading scalar of J=t^3*P*K^8*kappa on that locus.
    if kappa_numerator == 0:
        raise AssertionError("zero kappa numerator after determinant removal")

    return {
        "base": base,
        "residuals": residuals,
        "residual_metadata": residual_metadata,
        "determinant": D,
        "determinant_ratio": str(ratio),
        "lower_solution": {
            "m2": str(m2),
            "m1": str(m1),
            "m0": str(m0),
            "r1": str(r1),
            "r0": str(r0),
        },
        "kappa_numerator": kappa_numerator,
        "kappa_denominator": kappa_denominator,
        "kappa_removed_determinant_power": kappa_removed_power,
        "kappa_identity": (
            "coeff_t16((t-1)^2*N^3*K^2-T^2)=kappa*(-12)^8 "
            "on the differential locus"
        ),
    }


def polynomial_record(polynomial) -> dict[str, object]:
    return {
        "expression": str(polynomial),
        "total_degree": int(polynomial.total_degree()),
        "term_count": len(polynomial.dict()),
        "sha256": hashlib.sha256(str(polynomial).encode()).hexdigest(),
    }


def compute(output: Path) -> dict[str, object]:
    construction = construct_open_model_over_qq()
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
        raise AssertionError("a residual equation is not in J=(R1,R2)")

    resultant = R1b.resultant(R2b, r2b)
    plane = PolynomialRing(QQ, names=("k", "r3"), order="lex")
    kk, yy = plane.gens()
    resultant_plane = plane(str(resultant)).monic()
    factors = list(resultant_plane.factor())
    if len(factors) != 1 or int(factors[0][1]) != 1:
        raise AssertionError("the component plane resultant is not irreducible over QQ")

    ring = PolynomialRing(
        QQ,
        names=("u", "h", "r2", "r3", "k"),
        order="degrevlex",
    )
    u, h, r2, r3, k = ring.gens()
    embedding = base.hom([r2, r3, k], ring)
    open_product = (
        embedding(construction["determinant"])
        * embedding(construction["kappa_numerator"])
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
        test_ideal = ring.ideal(residuals + [saturation, u*target - 1])
        started = time.time()
        basis_raw = slimgb(test_ideal)
        elapsed = time.time() - started
        basis = [ring(value) for value in basis_raw]
        unit = any(
            value and int(value.total_degree()) == 0
            for value in basis
        )
        basis_text = "\n".join(str(value) for value in basis) + "\n"
        (output.parent / f"radical-membership-{name}.txt").write_text(
            basis_text,
            encoding="utf-8",
        )
        if not unit:
            raise AssertionError(("Rabinowitsch unit-ideal test failed", name))
        target_results[name] = {
            "rabinowitsch_unit_ideal": True,
            "elapsed_seconds": elapsed,
            "basis_count": len(basis),
            "basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
        }

    result: dict[str, object] = {
        "schema_version": 2,
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
        "constructor": {
            "base_field": "QQ",
            "method": "direct symbolic logarithmic-derivative construction",
            "residual_count": len(residuals),
            "residual_metadata": construction["residual_metadata"],
            "determinant": polynomial_record(construction["determinant"]),
            "determinant_ratio": construction["determinant_ratio"],
            "kappa_numerator": polynomial_record(
                construction["kappa_numerator"]
            ),
            "kappa_denominator": str(construction["kappa_denominator"]),
            "kappa_removed_determinant_power": construction[
                "kappa_removed_determinant_power"
            ],
            "kappa_identity": construction["kappa_identity"],
        },
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
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = compute(arguments.output)
    print(json.dumps({
        "constructor_base_field": result["constructor"]["base_field"],
        "radical_membership": result["radical_membership"],
        "set_theoretic_consequence": result["set_theoretic_consequence"],
        "record_sha256": result["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
