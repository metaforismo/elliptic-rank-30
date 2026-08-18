#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def primitive_integer_polynomial(expression, *generators):
    poly = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    _, primitive = poly.clear_denoms(convert=True)
    content, primitive = primitive.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def rational_function_data(expression, parameter):
    numerator, denominator = map(
        sp.factor, sp.cancel(expression).as_numer_denom()
    )
    return {
        "expression": str(sp.factor(expression)),
        "numerator_degree": int(sp.degree(numerator, parameter)),
        "denominator_degree": int(sp.degree(denominator, parameter)),
        "numerator_factorization": str(sp.factor(numerator)),
        "denominator_factorization": str(sp.factor(denominator)),
    }


def build_certificate():
    w = sp.symbols("w", real=True)
    x, y, z = sp.symbols("x y z")

    F = (
        81*x**5 + 1053*x**4*y - 1098*x**4
        + 13248*x**3*y - 11048*x**3
        + 9000*x**2*y**2 + 116688*x**2*y + 115872*x**2
        + 12000*x*y**2 + 87168*x*y + 112512*x
        + 40000*y**3 + 384000*y**2 + 979968*y + 813056
    )
    R1 = (
        -8496*x**2*y + 4374*x**2*z - 938*x**2
        - 2025*x*y**2 + 20552*x*y + 34992*x*z + 23236*x
        - 16425*y**3 + 29950*y**2 + 72900*y*z - 50456*y
        - 54216*z - 131648
    )
    R2 = (
        2832*x**3 + 18954*x**2*z - 65110*x**2
        - 8775*x*y**2 + 210520*x*y + 151632*x*z + 237884*x
        - 71175*y**3 - 169150*y**2 + 315900*y*z + 1062680*y
        + 1464264*z + 1262144
    )

    resultant_z = primitive_integer_polynomial(
        sp.resultant(R1, R2, z), x, y
    )
    F_primitive = primitive_integer_polynomial(F, x, y)
    if resultant_z != F_primitive:
        raise AssertionError("Res_z(R1,R2) did not equal the quintic")
    factorization = sp.factor_list(F)
    if len(factorization[1]) != 1 or factorization[1][0][1] != 1:
        raise AssertionError("the plane quintic is not irreducible over Q")

    k_w = sp.factor(
        4*(5*w**3 + 27*w**2 + 135*w + 81)
        / (3*(5*w - 9)*(w**2 + 3))
    )
    r3_w = sp.factor(
        -4*(w - 3)*(5*w**4 - 12*w**3 + 6*w**2 - 144*w + 81)
        / (3*(5*w - 9)*(w**2 + 3)**2)
    )
    r2_w = sp.factor(
        2*(w - 3)*(
            25*w**8 - 294*w**7 + 330*w**6 - 2790*w**5
            + 12420*w**4 - 162*w**3 + 58806*w**2
            - 24786*w + 19683
        )
        / (27*(5*w - 9)*(w**2 + 3)**4)
    )
    if sp.cancel(F.subs({x: k_w, y: r3_w})) != 0:
        raise AssertionError("the parametrization missed F")
    if sp.cancel(R1.subs({x: k_w, y: r3_w, z: r2_w})) != 0:
        raise AssertionError("the parametrization missed R1")
    if sp.cancel(R2.subs({x: k_w, y: r3_w, z: r2_w})) != 0:
        raise AssertionError("the parametrization missed R2")

    # Birational inverse on the affine plane component.
    x_equation = sp.together(x - k_w).as_numer_denom()[0]
    y_equation = sp.together(y - r3_w).as_numer_denom()[0]
    subresultants = sp.subresultants(x_equation, y_equation, w)
    inverse_linear = sp.Poly(sp.expand(subresultants[-2]), w)
    inverse_constant = primitive_integer_polynomial(
        subresultants[-1], x, y
    )
    if inverse_constant != F_primitive:
        raise AssertionError("the constant subresultant did not recover F")
    raw_A = sp.factor(inverse_linear.coeff_monomial(w))
    raw_B = sp.factor(inverse_linear.coeff_monomial(1))
    inverse_formula = sp.factor(-raw_B/raw_A)

    resultant_inverse_y = sp.factor(sp.resultant(F, raw_A, x))
    resultant_inverse_x = sp.factor(sp.resultant(F, raw_A, y))
    expected_y_support = (y + 4)**8 * (3*y + 4)**6
    expected_x_support = (x + 4)**8 * (3*x - 4)**6
    if sp.cancel(resultant_inverse_y/expected_y_support).free_symbols:
        raise AssertionError("unexpected affine inverse exceptional y-values")
    if sp.cancel(resultant_inverse_x/expected_x_support).free_symbols:
        raise AssertionError("unexpected affine inverse exceptional x-values")
    exceptional_affine = [
        {"x": "-4", "y": "-4", "z": "6", "parameter": "w=0"},
        {
            "x": "4/3",
            "y": "-4/3",
            "z": "10/27",
            "parameter": "w=infinity",
        },
    ]
    for point in exceptional_affine:
        px = sp.Rational(point["x"])
        py = sp.Rational(point["y"])
        pz = sp.Rational(point["z"])
        if F.subs({x: px, y: py}) != 0:
            raise AssertionError(("exceptional point is not on F", point))
        if R1.subs({x: px, y: py, z: pz}) != 0:
            raise AssertionError(("exceptional point is not on R1", point))
        if R2.subs({x: px, y: py, z: pz}) != 0:
            raise AssertionError(("exceptional point is not on R2", point))

    # The z-coordinate is rational on all rational parameter values in the
    # affine chart: the coefficient of z in R1 has no rational zero along P1.
    z_coefficient = sp.Poly(R1, z).coeff_monomial(z)
    z_coefficient_w = sp.factor(
        z_coefficient.subs({x: k_w, y: r3_w})
    )
    z_coefficient_numerator = sp.Poly(
        sp.together(z_coefficient_w).as_numer_denom()[0],
        w,
        domain=sp.QQ,
    )
    if z_coefficient_numerator.ground_roots():
        raise AssertionError(
            "the R1 z-coefficient has a rational parameter zero"
        )

    # Exact lower recurrence over QQ(w).
    field = sp.QQ.frac_field(w)
    zero = field.zero
    one = field.one

    def fe(expression):
        return field.from_sympy(sp.cancel(expression))

    kw = fe(k_w)
    r2w = fe(r2_w)
    r3w = fe(r3_w)

    def trim(values):
        result = list(values)
        while result and result[-1] == 0:
            result.pop()
        return result

    def add(left, right):
        size = max(len(left), len(right))
        return trim([
            (left[index] if index < len(left) else zero)
            + (right[index] if index < len(right) else zero)
            for index in range(size)
        ])

    def scale(values, scalar):
        scalar = scalar if hasattr(scalar, "parent") else fe(scalar)
        return trim([scalar*value for value in values])

    def subtract(left, right):
        return add(left, scale(right, -1))

    def multiply(left, right):
        if not left or not right:
            return []
        result = [zero]*(len(left) + len(right) - 1)
        for i, left_value in enumerate(left):
            for j, right_value in enumerate(right):
                result[i+j] += left_value*right_value
        return trim(result)

    def power(values, exponent):
        result = [one]
        base = list(values)
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = multiply(result, base)
            remaining //= 2
            if remaining:
                base = multiply(base, base)
        return result

    def derivative(values):
        return trim([
            fe(index)*values[index]
            for index in range(1, len(values))
        ])

    def divide_exact(dividend, divisor):
        numerator = trim(dividend)
        denominator = trim(divisor)
        if not denominator:
            raise ZeroDivisionError("zero polynomial")
        if len(numerator) < len(denominator):
            return [] if not numerator else None
        quotient = [zero]*(len(numerator) - len(denominator) + 1)
        inverse_lead = one/denominator[-1]
        while numerator and len(numerator) >= len(denominator):
            shift = len(numerator) - len(denominator)
            coefficient = numerator[-1]*inverse_lead
            quotient[shift] = coefficient
            for index, value in enumerate(denominator):
                numerator[index+shift] -= coefficient*value
            numerator = trim(numerator)
        return trim(quotient) if not numerator else None

    def coefficient(values, degree):
        return values[degree] if degree < len(values) else zero

    m2 = fe(8)*(-kw + fe(2)*r3w + fe(4))
    m1 = -(
        kw**2 + fe(14)*kw*r3w + fe(4)*kw + fe(24)*r2w
        + fe(13)*r3w**2 - fe(56)*r3w - fe(80)
    )/fe(3)

    def qbar(*, r0=zero, r1=zero, m0=zero):
        s = [fe(-1), one]
        R = [r0, r1, r2w, r3w, one]
        P = [zero] + R
        Q = add(derivative(P), scale(R, 3))
        K = [kw, fe(-12)]
        Kprime = [fe(-12)]
        M = [m0, m1, m2, fe(80)]
        N = add(power(Q, 2), multiply(P, M))
        first = multiply(
            multiply(N, K),
            subtract(multiply(s, Q), scale(P, 2)),
        )
        correction = scale(
            multiply(
                multiply(P, s),
                subtract(
                    multiply(derivative(N), K),
                    scale(multiply(N, Kprime), 2),
                ),
            ),
            3,
        )
        T = subtract(first, correction)
        E = subtract(
            add(
                subtract(
                    multiply(multiply(Q, T), K),
                    scale(
                        multiply(multiply(P, derivative(T)), K),
                        2,
                    ),
                ),
                scale(multiply(multiply(P, T), Kprime), 8),
            ),
            multiply(multiply(s, power(N, 2)), power(K, 2)),
        )
        quotient = divide_exact(E, P)
        if quotient is None:
            raise AssertionError("E lost its forced P factor")
        return quotient, {
            "R": R,
            "P": P,
            "Q": Q,
            "K": K,
            "M": M,
            "N": N,
            "T": T,
        }

    base_quotient, _ = qbar()
    degrees = [11, 10, 9]
    base_values = [
        coefficient(base_quotient, degree) for degree in degrees
    ]
    columns = []
    for name in ("m0", "r1", "r0"):
        trial, _ = qbar(**{name: one})
        columns.append([
            coefficient(trial, degree) - base_values[index]
            for index, degree in enumerate(degrees)
        ])
    matrix = [
        [columns[column][row] for column in range(3)]
        for row in range(3)
    ]
    rhs = [-value for value in base_values]

    def solve_three_by_three(A, b):
        work = [list(row) + [b[index]] for index, row in enumerate(A)]
        for column in range(3):
            pivot = next(
                row
                for row in range(column, 3)
                if work[row][column] != 0
            )
            if pivot != column:
                work[column], work[pivot] = work[pivot], work[column]
            inverse = one/work[column][column]
            work[column] = [value*inverse for value in work[column]]
            for row in range(3):
                if row == column or work[row][column] == 0:
                    continue
                factor = work[row][column]
                work[row] = [
                    work[row][index] - factor*work[column][index]
                    for index in range(4)
                ]
        return [work[index][-1] for index in range(3)]

    m0w, r1w, r0w = solve_three_by_three(matrix, rhs)
    final_quotient, pieces = qbar(r0=r0w, r1=r1w, m0=m0w)
    if any(final_quotient):
        raise AssertionError(
            "the parametrized differential equation is nonzero"
        )

    K2 = power(pieces["K"], 2)
    K4 = power(pieces["K"], 4)
    A = divide_exact(pieces["N"], K2)
    B = divide_exact(pieces["T"], K4)
    if A is None or B is None or len(A) != 7 or len(B) != 11:
        raise AssertionError("A/B reconstruction failed")
    if A[-1] != one or B[-1] != one:
        raise AssertionError("A/B lost monicity")

    p0_derived = sp.factor(field.to_sympy(-fe(4)*r0w/kw))
    p0_formula = sp.factor(
        -(w + 3)**11
        / (sp.Integer(243)*(w - 3)*(w**2 + 3)**5)
    )
    if sp.cancel(p0_derived - p0_formula) != 0:
        raise AssertionError("the p0 endpoint formula failed")

    def evaluate_at_one(values):
        result = zero
        for value in values:
            result += value
        return result

    B1_derived = sp.factor(field.to_sympy(
        -fe(2)
        * evaluate_at_one(pieces["N"])
        * evaluate_at_one(pieces["P"])
        / evaluate_at_one(pieces["K"])**3
    ))
    B1_formula = sp.factor(
        sp.Integer(2)**14*w**5*(w + 3)**7*(w**2 + 9)**10
        / (
            sp.Integer(3)**12*(w - 3)**2*(w**2 + 3)**15
        )
    )
    if sp.cancel(B1_derived - B1_formula) != 0:
        raise AssertionError("the B(1) endpoint formula failed")

    H = subtract(
        multiply(power([fe(-1), one], 2), power(A, 3)),
        power(B, 2),
    )
    target = [zero, zero, zero] + pieces["P"]
    if len(H) != len(target):
        raise AssertionError("unexpected H degree")
    kappa = H[-1]
    if H != scale(target, kappa):
        raise AssertionError("H is not kappa*t^3*P")
    kappa_derived = sp.factor(field.to_sympy(kappa))
    kappa_formula = sp.factor(
        sp.Integer(2)**22
        * w**8
        * (w + 3)**11
        * (5*w - 9)
        * (w**2 + 9)**16
        / (
            sp.Integer(3)**21
            * (w - 3)**4
            * (w**2 + 3)**24
        )
    )
    if sp.cancel(kappa_derived - kappa_formula) != 0:
        raise AssertionError("the discriminant scalar formula failed")

    determinant = (
        -2*x**2 - 37*x*y + 52*x + 108*z
        - 116*y**2 + 148*y - 320
    )
    determinant_w = sp.factor(
        determinant.subs({x: k_w, y: r3_w, z: r2_w})
    )
    determinant_numerator = sp.Poly(
        sp.together(determinant_w).as_numer_denom()[0],
        w,
        domain=sp.QQ,
    )
    if determinant_numerator.ground_roots() != {sp.Rational(3): 1}:
        raise AssertionError(
            "unexpected rational determinant-zero parameter"
        )

    i4_square_factor = sp.factor(
        (w + 3)**5/(9*(w**2 + 3)**2)
    )
    i4_square_class = sp.factor(
        (w + 3)/((w - 3)*(w**2 + 3))
    )
    iv_square_factor = sp.factor(
        sp.Integer(2)**7
        * w**2
        * (w + 3)**3
        * (w**2 + 9)**5
        / (
            sp.Integer(3)**6
            * (w - 3)
            * (w**2 + 3)**7
        )
    )
    iv_square_class = sp.factor(
        -2*w*(w + 3)/(w**2 + 3)
    )
    if sp.cancel(
        -3*p0_formula - i4_square_factor**2*i4_square_class
    ) != 0:
        raise AssertionError("I4 square-class factorization failed")
    if sp.cancel(
        -2*B1_formula - iv_square_factor**2*iv_square_class
    ) != 0:
        raise AssertionError("IV square-class factorization failed")

    i4_positive = sp.solve_univariate_inequality(
        i4_square_class > 0, w, relational=False
    )
    iv_positive = sp.solve_univariate_inequality(
        iv_square_class > 0, w, relational=False
    )
    positive_intersection = sp.Intersection(
        i4_positive, iv_positive
    )
    if positive_intersection != sp.EmptySet:
        raise AssertionError(
            "the two nonzero split square classes overlap"
        )

    infinity_point = {
        "k": str(sp.limit(k_w, w, sp.oo)),
        "r3": str(sp.limit(r3_w, w, sp.oo)),
        "r2": str(sp.limit(r2_w, w, sp.oo)),
        "p0": str(sp.limit(p0_formula, w, sp.oo)),
        "B_at_one": str(sp.limit(B1_formula, w, sp.oo)),
        "i4_target": str(sp.limit(-3*p0_formula, w, sp.oo)),
        "iv_target": str(sp.limit(-2*B1_formula, w, sp.oo)),
    }
    if sp.Rational(infinity_point["iv_target"]) >= 0:
        raise AssertionError(
            "the affine inverse exceptional infinity parameter "
            "has nonnegative IV target"
        )

    projective_variable = sp.symbols("Z")
    projective_leading = sp.factor(
        sp.Poly(F, x, y)
        .homogenize(projective_variable)
        .as_expr()
        .subs(projective_variable, 0)
    )
    if projective_leading != 81*x**4*(x + 13*y):
        raise AssertionError("unexpected projective points at infinity")

    payload = {
        "schema_version": 1,
        "certificate_id": (
            "rank17_iv_rational_component_split_obstruction"
        ),
        "truth_status": (
            "CERTIFIED exact rational parametrization and real-sign "
            "obstruction for one irreducible additive-IV surface "
            "component; completeness of this component for the entire "
            "saturated IV surface locus remains a separate elimination "
            "gate; no height-79/12 section or rank-30 conclusion"
        ),
        "solved_a": False,
        "solved_b": False,
        "unconditional_global_rank_lower_bound": 29,
        "plane_component": {
            "quintic": str(F_primitive.as_expr()),
            "irreducible_over_Q": True,
            "complete_intersection_relations": [str(R1), str(R2)],
            "resultant_z_equals_quintic": True,
            "rational_parametrization": {
                "k": str(k_w),
                "r3": str(r3_w),
                "r2": str(r2_w),
            },
            "parametrization_identities_verified": [
                "F=0",
                "R1=0",
                "R2=0",
            ],
            "inverse_linear_subresultant": {
                "formula": f"w={inverse_formula}",
                "coefficient_of_w": str(raw_A),
                "constant_term": str(raw_B),
                "affine_exceptional_points": exceptional_affine,
                "exceptional_support_resultant_in_y": str(
                    resultant_inverse_y
                ),
                "exceptional_support_resultant_in_x": str(
                    resultant_inverse_x
                ),
            },
            "projective_infinity_factor": str(projective_leading),
            "z_lift_coefficient_along_parameter": str(
                z_coefficient_w
            ),
            "z_lift_coefficient_has_no_rational_zero": True,
        },
        "logarithmic_derivative_reconstruction": {
            "differential_equation_zero_over_Qw": True,
            "r0": rational_function_data(field.to_sympy(r0w), w),
            "r1": rational_function_data(field.to_sympy(r1w), w),
            "m0": rational_function_data(field.to_sympy(m0w), w),
            "p0": rational_function_data(p0_formula, w),
            "B_at_one": rational_function_data(B1_formula, w),
            "discriminant_scalar_kappa": rational_function_data(
                kappa_formula, w
            ),
            "lower_determinant_along_parameter": str(determinant_w),
            "only_rational_determinant_zero_parameter": "w=3",
        },
        "split_square_classes": {
            "I4_target": "-3*p0",
            "I4_exact_factorization": (
                f"({i4_square_factor})^2 * ({i4_square_class})"
            ),
            "I4_nonzero_positive_set": str(i4_positive),
            "IV_target": "-2*B(1)",
            "IV_exact_factorization": (
                f"({iv_square_factor})^2 * ({iv_square_class})"
            ),
            "IV_nonzero_positive_set": str(iv_positive),
            "positive_set_intersection": str(positive_intersection),
            "real_sign_consequence": (
                "No finite nondegenerate rational parameter can make "
                "both split targets rational squares."
            ),
        },
        "exceptional_parameter_checks": {
            "w=-3": "p0=B(1)=0; degenerate fibre data",
            "w=0": (
                "kappa=B(1)=0; zero-discriminant boundary"
            ),
            "w=3": (
                "lower determinant vanishes; p0 and B(1) have "
                "poles, so this normalization reaches only a "
                "projective boundary point"
            ),
            "w=9/5": (
                "affine moduli coordinates have a projective pole "
                "and kappa=0"
            ),
            "w=infinity": infinity_point,
        },
        "mathematical_consequence": (
            "The nondegenerate Q-rational points of this rational "
            "component do not realize both required split fibres.  "
            "Hence this component cannot supply the sought split "
            "I12+I4+IV rank-17 seed."
        ),
        "limitations": [
            (
                "The determinant-zero chart remains separate.  Its "
                "only rational parameter on this normalization is "
                "w=3, where the original affine surface coordinates "
                "have poles."
            ),
            (
                "The certificate does not yet prove that this quintic "
                "is the only irreducible component of the fully "
                "saturated characteristic-zero IV surface locus."
            ),
            (
                "The section incidence and height-79/12 conditions "
                "are not imposed because the split-fibre obstruction "
                "occurs earlier."
            ),
            "The semistable I3 realization remains a separate branch.",
            "The unconditional global rank lower bound remains 29.",
        ],
    }
    payload["record_sha256"] = canonical_hash(payload)
    return payload


def markdown(payload):
    return "\n".join([
        "# Rational additive-IV component: split-fibre obstruction",
        "",
        "```text",
        "SOLVED-A: false",
        "SOLVED-B: false",
        "unconditional global lower bound: rank E(Q) >= 29",
        "```",
        "",
        payload["truth_status"],
        "",
        "## Certified consequence",
        "",
        payload["mathematical_consequence"],
        "",
        "## Exact square classes",
        "",
        (
            "- I4 positivity: `"
            f"{payload['split_square_classes']['I4_nonzero_positive_set']}`"
        ),
        (
            "- IV positivity: `"
            f"{payload['split_square_classes']['IV_nonzero_positive_set']}`"
        ),
        (
            "- intersection: `"
            f"{payload['split_square_classes']['positive_set_intersection']}`"
        ),
        "",
        "## Boundary",
        "",
        *[f"- {item}" for item in payload["limitations"]],
        "",
        f"Record SHA-256: `{payload['record_sha256']}`",
        "",
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    payload = build_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if arguments.markdown_output:
        arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.markdown_output.write_text(
            markdown(payload), encoding="utf-8"
        )
    if arguments.compare:
        committed = json.loads(
            arguments.compare.read_text(encoding="utf-8")
        )
        if committed != payload:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
    print(json.dumps({
        "certificate_id": payload["certificate_id"],
        "positive_set_intersection": payload[
            "split_square_classes"
        ]["positive_set_intersection"],
        "record_sha256": payload["record_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
