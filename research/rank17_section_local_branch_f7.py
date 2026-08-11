#!/usr/bin/env python3
"""Exact local F_7 certificate for the rank-17 section-incidence point.

The input point is the unique target-component section found on split surface 3
in the exhaustive all-chart search.  This script works only over F_7 and proves:

* the 19 surface-plus-section equations vanish at the point;
* their 19x17 Jacobian has rank 14;
* after solving 14 pivot equations formally, the residual tangent cone in the
  three free coordinates a=y1-1, u=y2-4, v=y3-2 is

      4 (u+v)^4,
      a^3 + a^2(u+v) - 2a(u+v)^2 - (u+v)^3;

* the reduced tangent direction is a=0, v=-u;
* among quadratic corrections a=A2*tau^2, v=-tau+V2*tau^2, the unique pair
  satisfying all residual equations through order 8 is (A2,V2)=(2,5);
* that second-order arc satisfies every residual through order 10, but with no
  further free coefficients the first obstruction is coefficient 6 at order 12
  in residual equation 15.

This is a certified local computation in characteristic 7.  It is not a
characteristic-zero section, a p-adic lift, or a rank-30 certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable

P = 7
MAX_DEGREE = 4
VARIABLES = [
    "p0", "p1", "p2", "p3", "q0", "q1", "q2", "s",
    "x0", "x1", "x2", "x3", "x4", "y0", "y1", "y2", "y3",
]
SEED = [2, 2, 1, 0, 2, 2, 2, 0, 5, 6, 2, 0, 6, 0, 1, 4, 2]
PIVOT_COLUMNS = list(range(14))
PIVOT_ROWS = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14]
RESIDUAL_ROWS = [6, 15, 16, 17, 18]


def mod(value: int) -> int:
    return value % P


class Dual:
    __slots__ = ("value", "gradient")

    def __init__(self, value: int, gradient: Iterable[int]):
        self.value = mod(value)
        self.gradient = tuple(mod(entry) for entry in gradient)

    def _coerce(self, other: object) -> "Dual":
        if isinstance(other, Dual):
            return other
        return Dual(int(other), [0] * len(self.gradient))

    def __add__(self, other: object) -> "Dual":
        rhs = self._coerce(other)
        return Dual(
            self.value + rhs.value,
            [a + b for a, b in zip(self.gradient, rhs.gradient, strict=True)],
        )

    __radd__ = __add__

    def __neg__(self) -> "Dual":
        return Dual(-self.value, [-entry for entry in self.gradient])

    def __sub__(self, other: object) -> "Dual":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "Dual":
        return self._coerce(other) - self

    def __mul__(self, other: object) -> "Dual":
        rhs = self._coerce(other)
        return Dual(
            self.value * rhs.value,
            [
                self.value * db + rhs.value * da
                for da, db in zip(self.gradient, rhs.gradient, strict=True)
            ],
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "Dual":
        rhs = self._coerce(other)
        if any(rhs.gradient):
            raise ValueError("dual division is used only by constants here")
        inverse = pow(rhs.value, -1, P)
        return Dual(
            self.value * inverse,
            [entry * inverse for entry in self.gradient],
        )

    def __pow__(self, exponent: int) -> "Dual":
        if exponent < 0:
            raise ValueError("negative powers are not used")
        result = Dual(1, [0] * len(self.gradient))
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result


class Series:
    """Sparse truncated multivariate series over F_7."""

    __slots__ = ("coefficients", "max_degree", "variable_count")

    def __init__(
        self,
        value: int | dict[tuple[int, ...], int],
        max_degree: int,
        variable_count: int,
    ):
        self.max_degree = max_degree
        self.variable_count = variable_count
        if isinstance(value, dict):
            self.coefficients = {
                tuple(monomial): mod(coefficient)
                for monomial, coefficient in value.items()
                if coefficient % P and sum(monomial) <= max_degree
            }
        else:
            self.coefficients = (
                {(0,) * variable_count: mod(value)} if value % P else {}
            )

    def _coerce(self, other: object) -> "Series":
        if isinstance(other, Series):
            if (
                other.max_degree != self.max_degree
                or other.variable_count != self.variable_count
            ):
                raise ValueError("incompatible formal-series rings")
            return other
        return Series(int(other), self.max_degree, self.variable_count)

    def copy(self) -> "Series":
        return Series(dict(self.coefficients), self.max_degree, self.variable_count)

    def __add__(self, other: object) -> "Series":
        rhs = self._coerce(other)
        result = dict(self.coefficients)
        for monomial, coefficient in rhs.coefficients.items():
            new_value = mod(result.get(monomial, 0) + coefficient)
            if new_value:
                result[monomial] = new_value
            else:
                result.pop(monomial, None)
        return Series(result, self.max_degree, self.variable_count)

    __radd__ = __add__

    def __neg__(self) -> "Series":
        return Series(
            {monomial: -coefficient for monomial, coefficient in self.coefficients.items()},
            self.max_degree,
            self.variable_count,
        )

    def __sub__(self, other: object) -> "Series":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "Series":
        return self._coerce(other) - self

    def __mul__(self, other: object) -> "Series":
        rhs = self._coerce(other)
        result: dict[tuple[int, ...], int] = {}
        for left_monomial, left_coefficient in self.coefficients.items():
            for right_monomial, right_coefficient in rhs.coefficients.items():
                monomial = tuple(
                    a + b
                    for a, b in zip(
                        left_monomial, right_monomial, strict=True
                    )
                )
                if sum(monomial) > self.max_degree:
                    continue
                result[monomial] = mod(
                    result.get(monomial, 0)
                    + left_coefficient * right_coefficient
                )
        return Series(result, self.max_degree, self.variable_count)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "Series":
        rhs = self._coerce(other)
        zero = (0,) * self.variable_count
        if set(rhs.coefficients) != {zero}:
            raise ValueError("only constant formal-series division is used")
        inverse = pow(rhs.coefficients[zero], -1, P)
        return Series(
            {
                monomial: coefficient * inverse
                for monomial, coefficient in self.coefficients.items()
            },
            self.max_degree,
            self.variable_count,
        )

    def __pow__(self, exponent: int) -> "Series":
        if exponent < 0:
            raise ValueError("negative powers are not used")
        result = Series(1, self.max_degree, self.variable_count)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def coefficient(self, monomial: tuple[int, ...]) -> int:
        return self.coefficients.get(tuple(monomial), 0)

    def homogeneous(self, degree: int) -> dict[tuple[int, ...], int]:
        return {
            monomial: coefficient
            for monomial, coefficient in self.coefficients.items()
            if sum(monomial) == degree
        }


def convolution(values, other, length):
    zero = values[0] * 0
    result = [zero for _ in range(length)]
    for left_index, left_value in enumerate(values):
        for right_index, right_value in enumerate(other):
            if left_index + right_index < length:
                result[left_index + right_index] = (
                    result[left_index + right_index]
                    + left_value * right_value
                )
    return result


def equations(values):
    p0, p1, p2, p3, q0, q1, q2, s, *rest = values
    zero = p0 * 0
    one = zero + 1
    x_values = rest[:5] + [zero] * 4
    y_values = rest[5:9] + [zero] * 9

    c4 = [zero] * 9
    c4[0] = p0 * p0
    c4[1] = 2 * p0 * p1
    c4[2] = 2 * p0 * p2 + p1 * p1
    c4[3] = 2 * p0 * p3 + 2 * p1 * p2
    c4[4] = (
        -3 - 15 * p0 * p0 - 20 * p0 * p1 - 12 * p0 * p2
        - 6 * p0 * p3 - 6 * p1 * p1 - 6 * p1 * p2
        + 15 * q0 * q0 - 10 * q0 * q1 + 2 * q0 * q2
        + q1 * q1 - s
    )
    c4[5] = (
        8 + 24 * p0 * p0 + 30 * p0 * p1 + 16 * p0 * p2
        + 6 * p0 * p3 + 8 * p1 * p1 + 6 * p1 * p2
        - 24 * q0 * q0 + 18 * q0 * q1 - 4 * q0 * q2
        - 2 * q1 * q1 + 3 * s
    )
    c4[6] = (
        -6 - 10 * p0 * p0 - 12 * p0 * p1 - 6 * p0 * p2
        - 2 * p0 * p3 - 3 * p1 * p1 - 2 * p1 * p2
        + 10 * q0 * q0 - 8 * q0 * q1 + 2 * q0 * q2
        + q1 * q1 - 3 * s
    )
    c4[7] = s
    c4[8] = one

    reversed_c4 = list(reversed(c4))
    c4_cubed = convolution(
        convolution(reversed_c4, reversed_c4, 12),
        reversed_c4,
        12,
    )
    square_root = [one]
    for order in range(1, 12):
        correction = zero
        for index in range(1, order):
            correction = (
                correction
                + square_root[index] * square_root[order - index]
            )
        square_root.append((c4_cubed[order] - correction) / 2)

    c6 = [zero] * 13
    c6[0] = p0**3
    c6[12] = one
    for order in range(1, 12):
        c6[12 - order] = square_root[order]

    result = [
        c6[1] - 3 * p0 * p0 * p1,
        c6[2] - 3 * (p0 * p0 * p2 + p0 * p1 * p1),
        c6[3]
        - (3 * p0 * p0 * p3 + 6 * p0 * p1 * p2 + p1**3),
        sum(c6, zero) - q0**3,
        sum((index * c6[index] for index in range(13)), zero)
        - 3 * q0 * q0 * q1,
        sum(
            (
                (index * (index - 1) // 2) * c6[index]
                for index in range(13)
            ),
            zero,
        )
        - 3 * (q0 * q0 * q2 + q0 * q1 * q1),
    ]

    y_squared = convolution(y_values, y_values, 13)
    x_squared = convolution(x_values, x_values, 13)
    x_cubed = convolution(x_squared, x_values, 13)
    c4_times_x = convolution(c4, x_values, 13)
    result.extend(
        y_squared[index]
        - x_cubed[index]
        + 3 * c4_times_x[index]
        + 2 * c6[index]
        for index in range(13)
    )
    return result


def rref(matrix):
    work = [[mod(value) for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    pivot_columns = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (
                row
                for row in range(pivot_row, row_count)
                if work[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, P)
        work[pivot_row] = [
            mod(value * inverse) for value in work[pivot_row]
        ]
        for row in range(row_count):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    mod(value - factor * pivot_value)
                    for value, pivot_value in zip(
                        work[row], work[pivot_row], strict=True
                    )
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return work, pivot_columns


def solve_square(matrix, right_hand_side):
    size = len(matrix)
    augmented = [
        [mod(value) for value in matrix[row]]
        + [mod(right_hand_side[row])]
        for row in range(size)
    ]
    for column in range(size):
        selected = next(
            row
            for row in range(column, size)
            if augmented[row][column]
        )
        augmented[column], augmented[selected] = (
            augmented[selected], augmented[column]
        )
        inverse = pow(augmented[column][column], -1, P)
        augmented[column] = [
            mod(value * inverse) for value in augmented[column]
        ]
        for row in range(size):
            if row != column and augmented[row][column]:
                factor = augmented[row][column]
                augmented[row] = [
                    mod(value - factor * pivot_value)
                    for value, pivot_value in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
    return [augmented[row][-1] for row in range(size)]


def monomials_of_degree(variable_count, degree):
    if variable_count == 1:
        return [(degree,)]
    result = []

    def recurse(prefix, remaining, positions):
        if positions == 1:
            result.append(tuple(prefix + [remaining]))
            return
        for exponent in range(remaining + 1):
            recurse(prefix + [exponent], remaining - exponent, positions - 1)

    recurse([], degree, variable_count)
    return result


def jacobian_at_seed():
    dimension = len(SEED)
    dual_values = []
    for index, value in enumerate(SEED):
        gradient = [0] * dimension
        gradient[index] = 1
        dual_values.append(Dual(value, gradient))
    evaluated = equations(dual_values)
    if any(item.value for item in evaluated):
        raise AssertionError("modular point is not on the incidence scheme")
    return [list(item.gradient) for item in evaluated]


def independent_rows(matrix, target_rank):
    selected = []
    current = []
    rank = 0
    for index, row in enumerate(matrix):
        candidate = current + [row]
        _, pivot_columns = rref(candidate)
        if len(pivot_columns) > rank:
            selected.append(index)
            current = candidate
            rank += 1
        if rank == target_rank:
            return selected
    raise AssertionError("not enough independent rows")


def formal_implicit_solution(jacobian):
    variables = [Series(value, MAX_DEGREE, 3) for value in SEED]
    free_monomials = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for column, monomial in zip([14, 15, 16], free_monomials, strict=True):
        variables[column].coefficients[monomial] = 1

    pivot_matrix = [
        [jacobian[row][column] for column in PIVOT_COLUMNS]
        for row in PIVOT_ROWS
    ]
    _, pivots = rref(pivot_matrix)
    if pivots != list(range(14)):
        raise AssertionError("selected pivot matrix is singular")

    for degree in range(1, MAX_DEGREE + 1):
        evaluated = equations(variables)
        for monomial in monomials_of_degree(3, degree):
            right_hand_side = [
                -evaluated[row].coefficient(monomial)
                for row in PIVOT_ROWS
            ]
            correction = solve_square(pivot_matrix, right_hand_side)
            for column, coefficient in zip(
                PIVOT_COLUMNS, correction, strict=True
            ):
                if coefficient:
                    current = variables[column].coefficient(monomial)
                    new_value = mod(current + coefficient)
                    if new_value:
                        variables[column].coefficients[monomial] = new_value
                    else:
                        variables[column].coefficients.pop(monomial, None)
        verified = equations(variables)
        for row in PIVOT_ROWS:
            if any(
                coefficient
                for monomial, coefficient in verified[row].coefficients.items()
                if sum(monomial) <= degree
            ):
                raise AssertionError(("pivot equation not solved", degree, row))
    return variables, equations(variables)


def expand_expected_cubic():
    # a^3 + a^2*w - 2*a*w^2 - w^3, w=u+v.
    result: dict[tuple[int, int, int], int] = {}

    def add_term(monomial, coefficient):
        result[monomial] = mod(result.get(monomial, 0) + coefficient)
        if not result[monomial]:
            result.pop(monomial)

    add_term((3, 0, 0), 1)
    add_term((2, 1, 0), 1)
    add_term((2, 0, 1), 1)
    add_term((1, 2, 0), -2)
    add_term((1, 1, 1), -4)
    add_term((1, 0, 2), -2)
    add_term((0, 3, 0), -1)
    add_term((0, 2, 1), -3)
    add_term((0, 1, 2), -3)
    add_term((0, 0, 3), -1)
    return result


def expand_expected_quartic():
    # 4*(u+v)^4.
    return {
        (0, 4, 0): 4,
        (0, 3, 1): mod(16),
        (0, 2, 2): mod(24),
        (0, 1, 3): mod(16),
        (0, 0, 4): 4,
    }


def substitute_arc(series, a2, v2, order):
    result = [0] * (order + 1)
    # a=a2*tau^2, u=tau, v=-tau+v2*tau^2.
    for (a_exp, u_exp, v_exp), coefficient in series.coefficients.items():
        a_poly = [0] * (order + 1)
        if a_exp == 0:
            a_poly[0] = 1
        elif 2 * a_exp <= order:
            a_poly[2 * a_exp] = pow(a2, a_exp, P)

        u_poly = [0] * (order + 1)
        if u_exp <= order:
            u_poly[u_exp] = 1

        v_base = [0] * (order + 1)
        if order >= 1:
            v_base[1] = -1 % P
        if order >= 2:
            v_base[2] = v2 % P

        def multiply_poly(left, right):
            product = [0] * (order + 1)
            for left_degree, left_value in enumerate(left):
                if not left_value:
                    continue
                for right_degree, right_value in enumerate(right):
                    if left_degree + right_degree > order:
                        break
                    product[left_degree + right_degree] = mod(
                        product[left_degree + right_degree]
                        + left_value * right_value
                    )
            return product

        def power_poly(base, exponent):
            value = [1] + [0] * order
            current = base
            remaining = exponent
            while remaining:
                if remaining & 1:
                    value = multiply_poly(value, current)
                current = multiply_poly(current, current)
                remaining //= 2
            return value

        term = multiply_poly(a_poly, u_poly)
        term = multiply_poly(term, power_poly(v_base, v_exp))
        for degree, value in enumerate(term):
            result[degree] = mod(result[degree] + coefficient * value)
    return result


def univariate_arc_obstruction(jacobian, order=12, a2=2, v2=5):
    # Re-solve the 14 pivot equations as true univariate series, so this check
    # does not rely on the total-degree-10 multivariate truncation.
    variables = [Series(value, order, 1) for value in SEED]
    tau = (1,)
    if a2:
        variables[14].coefficients[(2,)] = a2 % P
    variables[15].coefficients[tau] = 1
    variables[16].coefficients[tau] = -1 % P
    if v2:
        variables[16].coefficients[(2,)] = v2 % P

    pivot_matrix = [
        [jacobian[row][column] for column in PIVOT_COLUMNS]
        for row in PIVOT_ROWS
    ]
    for degree in range(1, order + 1):
        evaluated = equations(variables)
        monomial = (degree,)
        correction = solve_square(
            pivot_matrix,
            [-evaluated[row].coefficient(monomial) for row in PIVOT_ROWS],
        )
        for column, coefficient in zip(PIVOT_COLUMNS, correction, strict=True):
            if coefficient:
                current = variables[column].coefficient(monomial)
                new_value = mod(current + coefficient)
                if new_value:
                    variables[column].coefficients[monomial] = new_value
                else:
                    variables[column].coefficients.pop(monomial, None)
    evaluated = equations(variables)
    residuals = {
        str(row): [
            [degree, evaluated[row].coefficient((degree,))]
            for degree in range(order + 1)
            if evaluated[row].coefficient((degree,))
        ]
        for row in RESIDUAL_ROWS
    }
    return variables, residuals


def canonical_hash(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def compute_certificate():
    direct = equations(SEED)
    if any(item % P for item in direct):
        raise AssertionError("seed equations do not vanish")

    jacobian = jacobian_at_seed()
    _, pivot_columns = rref(jacobian)
    if pivot_columns != PIVOT_COLUMNS:
        raise AssertionError(("unexpected Jacobian pivots", pivot_columns))
    selected_rows = independent_rows(
        [[row[column] for column in PIVOT_COLUMNS] for row in jacobian],
        14,
    )
    if selected_rows != PIVOT_ROWS:
        raise AssertionError(("unexpected pivot rows", selected_rows))

    formal_variables, evaluated = formal_implicit_solution(jacobian)
    residual_leading = {}
    for row in RESIDUAL_ROWS:
        nonzero_degrees = sorted(
            {sum(monomial) for monomial in evaluated[row].coefficients}
        )
        degree = nonzero_degrees[0] if nonzero_degrees else None
        residual_leading[str(row)] = {
            "degree": degree,
            "coefficients": [
                {"monomial_a_u_v": list(monomial), "coefficient": coefficient}
                for monomial, coefficient in sorted(
                    evaluated[row].homogeneous(degree).items()
                )
            ] if degree is not None else [],
        }

    quartic = evaluated[6].homogeneous(4)
    cubic = evaluated[15].homogeneous(3)
    if quartic != expand_expected_quartic():
        raise AssertionError(("unexpected quartic tangent equation", quartic))
    if cubic != expand_expected_cubic():
        raise AssertionError(("unexpected cubic tangent equation", cubic))

    # v=-u kills the quartic and leaves exactly a^3 in the cubic.
    cubic_on_v_minus_u = [0] * 4
    for (a_exp, u_exp, v_exp), coefficient in cubic.items():
        total_u = u_exp + v_exp
        sign = -1 if v_exp % 2 else 1
        if total_u == 0:
            cubic_on_v_minus_u[a_exp] = mod(
                cubic_on_v_minus_u[a_exp] + sign * coefficient
            )
        elif a_exp == 0:
            # These pure-u terms must cancel coefficient-by-coefficient.
            cubic_on_v_minus_u[0] = mod(
                cubic_on_v_minus_u[0] + sign * coefficient
            )
        else:
            # Store a mixed term in a separate sentinel by direct expansion.
            pass
    # Direct exhaustive symbolic evaluation is simpler and fail-closed.
    for a_value in range(P):
        for u_value in range(P):
            cubic_value = sum(
                coefficient
                * pow(a_value, a_exp, P)
                * pow(u_value, u_exp, P)
                * pow(-u_value % P, v_exp, P)
                for (a_exp, u_exp, v_exp), coefficient in cubic.items()
            ) % P
            if cubic_value != pow(a_value, 3, P):
                raise AssertionError("cubic restriction is not a^3")
            quartic_value = sum(
                coefficient
                * pow(a_value, a_exp, P)
                * pow(u_value, u_exp, P)
                * pow(-u_value % P, v_exp, P)
                for (a_exp, u_exp, v_exp), coefficient in quartic.items()
            ) % P
            if quartic_value:
                raise AssertionError("quartic does not vanish on v=-u")

    second_order_solutions = []
    for a2 in range(P):
        for v2 in range(P):
            _variables, residuals = univariate_arc_obstruction(
                jacobian, 8, a2=a2, v2=v2
            )
            if all(not values for values in residuals.values()):
                second_order_solutions.append([a2, v2])
    if second_order_solutions != [[2, 5]]:
        raise AssertionError(
            ("unexpected quadratic corrections", second_order_solutions)
        )

    _variables10, residuals10 = univariate_arc_obstruction(
        jacobian, 10, a2=2, v2=5
    )
    if any(residuals10.values()):
        raise AssertionError(
            ("selected arc obstructed before order 11", residuals10)
        )

    univariate_variables, residual_obstructions = univariate_arc_obstruction(
        jacobian, 12, a2=2, v2=5
    )
    expected_obstructions = {
        "6": [],
        "15": [[12, 6]],
        "16": [],
        "17": [],
        "18": [],
    }
    if residual_obstructions != expected_obstructions:
        raise AssertionError(
            ("unexpected order-12 obstruction", residual_obstructions)
        )

    payload = {
        "schema_version": 1,
        "certificate_id": "rank17_section_local_branch_f7",
        "truth_status": (
            "CERTIFIED local F7 tangent-cone and second-order arc; "
            "full p-adic and characteristic-zero lift unresolved"
        ),
        "modular_incidence_point": {
            "prime": P,
            "variables": VARIABLES,
            "values": SEED,
            "surface_parameters": {
                name: SEED[index] for index, name in enumerate(VARIABLES[:8])
            },
            "section": {
                "D": [1, 0, 0],
                "X_coefficients_ascending": SEED[8:13] + [0, 0, 0, 0],
                "Y_coefficients_ascending": SEED[13:17] + [0] * 9,
                "required_component_tangent_ratios": [1, 1],
            },
            "all_19_equations_vanish": True,
        },
        "jacobian": {
            "shape": [19, 17],
            "rank_over_F7": 14,
            "matrix": jacobian,
            "pivot_columns": PIVOT_COLUMNS,
            "pivot_variables": [VARIABLES[index] for index in PIVOT_COLUMNS],
            "free_columns": [14, 15, 16],
            "free_coordinates": ["a=y1-1", "u=y2-4", "v=y3-2"],
            "independent_pivot_rows": PIVOT_ROWS,
        },
        "formal_elimination": {
            "truncation_total_degree": MAX_DEGREE,
            "pivot_equations_solved": PIVOT_ROWS,
            "residual_equations": RESIDUAL_ROWS,
            "leading_terms": residual_leading,
            "factored_tangent_equations": {
                "row_6_degree_4": "4*(u+v)^4",
                "row_15_degree_3": (
                    "a^3+a^2*(u+v)-2*a*(u+v)^2-(u+v)^3"
                ),
                "row_15_restricted_to_v=-u": "a^3",
            },
            "reduced_tangent_direction": {
                "a": 0,
                "u": 1,
                "v": -1,
                "interpretation": (
                    "The reduced tangent cone is the line a=0, u+v=0."
                ),
            },
        },
        "quadratic_arc": {
            "parameter": "tau",
            "free_coordinate_expansions": {
                "y1": "1+2*tau^2",
                "y2": "4+tau",
                "y3": "2-tau-2*tau^2 (mod 7)",
            },
            "quadratic_correction_candidates_tested": P * P,
            "unique_solution_A2_V2": [2, 5],
            "all_residuals_zero_through_order": 10,
            "true_univariate_reelimination_order": 12,
            "first_obstruction_without_higher_corrections": {
                "residual_row": 15,
                "tau_degree": 12,
                "coefficient_mod_7": 6,
            },
        },
        "limitations": [
            (
                "A reduced tangent line and a unique quadratic correction do "
                "not by themselves prove an infinite formal branch."
            ),
            (
                "The order-12 obstruction shows that higher corrections or a "
                "deflated local system are required."
            ),
            (
                "No p-adic point, rational section, rank-17 fibration, or "
                "rank-30 curve is certified by this artifact."
            ),
        ],
        "implementation": {
            "language": "Python standard library",
            "script": "research/rank17_section_local_branch_f7.py",
        },
    }
    payload["certificate_sha256"] = canonical_hash(payload)
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", type=Path)
    arguments = parser.parse_args()
    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text())
        if committed != certificate:
            raise AssertionError(
                f"certificate mismatch: {arguments.compare}"
            )
    print(
        json.dumps(
            {
                "truth_status": certificate["truth_status"],
                "jacobian_rank": certificate["jacobian"]["rank_over_F7"],
                "quadratic_correction": certificate["quadratic_arc"]["unique_solution_A2_V2"],
                "first_obstruction": certificate["quadratic_arc"]["first_obstruction_without_higher_corrections"],
                "certificate_sha256": certificate["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
