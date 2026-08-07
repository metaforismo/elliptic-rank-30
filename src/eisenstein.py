"""Exact arithmetic in Q(zeta_3) for cyclic cubic packet certificates."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence


@dataclass(frozen=True)
class Eisenstein:
    """The element a+b*zeta, where zeta^2+zeta+1=0."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __init__(self, a=0, b=0):
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))

    def __add__(self, other):
        other = coerce(other)
        return Eisenstein(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-coerce(other))

    def __rsub__(self, other):
        return coerce(other) - self

    def __mul__(self, other):
        other = coerce(other)
        return Eisenstein(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    __rmul__ = __mul__

    def conjugate(self):
        # zeta_bar=zeta^2=-1-zeta.
        return Eisenstein(self.a - self.b, -self.b)

    def norm(self):
        return self.a * self.a - self.a * self.b + self.b * self.b

    def inverse(self):
        if not self:
            raise ZeroDivisionError("division by zero in Q(zeta_3)")
        conjugate = self.conjugate()
        norm = self.norm()
        return Eisenstein(conjugate.a / norm, conjugate.b / norm)

    def __truediv__(self, other):
        return self * coerce(other).inverse()

    def __rtruediv__(self, other):
        return coerce(other) / self

    def __bool__(self):
        return bool(self.a or self.b)

    def is_integral(self):
        return self.a.denominator == 1 and self.b.denominator == 1

    def as_pair(self):
        return self.a, self.b

    def __str__(self):
        if not self.b:
            return str(self.a)
        if not self.a:
            return f"{self.b}*zeta"
        sign = "+" if self.b >= 0 else "-"
        return f"{self.a}{sign}{abs(self.b)}*zeta"


def coerce(value):
    return value if isinstance(value, Eisenstein) else Eisenstein(value)


def determinant(matrix: Sequence[Sequence[Eisenstein]]) -> Eisenstein:
    """Exact Gaussian determinant over Q(zeta_3)."""

    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("matrix must be square")
    work = [[coerce(entry) for entry in row] for row in matrix]
    result = Eisenstein(1)
    sign = 1
    for column in range(n):
        pivot = next((row for row in range(column, n) if work[row][column]), None)
        if pivot is None:
            return Eisenstein(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, n):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for j in range(column + 1, n):
                work[row][j] -= factor * work[column][j]
            work[row][column] = Eisenstein(0)
    return result if sign > 0 else -result


def is_hermitian(matrix: Sequence[Sequence[Eisenstein]]) -> bool:
    n = len(matrix)
    return all(
        coerce(matrix[i][j]) == coerce(matrix[j][i]).conjugate()
        for i in range(n)
        for j in range(n)
    )


def trace_form_gram(matrix: Sequence[Sequence[Eisenstein]]):
    """Return the Z-bilinear trace-form Gram matrix in bases (1,zeta)."""

    n = len(matrix)
    if not is_hermitian(matrix):
        raise ValueError("Hermitian matrix required")
    basis = (Eisenstein(1), Eisenstein(0, 1))
    gram = [[Fraction(0) for _ in range(2 * n)] for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            hij = coerce(matrix[i][j])
            for r, left in enumerate(basis):
                for s, right in enumerate(basis):
                    value = left.conjugate() * hij * right
                    # Field trace of a+b*zeta is 2a-b.
                    gram[2 * i + r][2 * j + s] = 2 * value.a - value.b
    return gram


def determinant_rational(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    n = len(matrix)
    work = [[Fraction(entry) for entry in row] for row in matrix]
    result = Fraction(1)
    sign = 1
    for column in range(n):
        pivot = next((row for row in range(column, n) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, n):
            factor = work[row][column] / pivot_value
            for j in range(column + 1, n):
                work[row][j] -= factor * work[column][j]
            work[row][column] = Fraction(0)
    return result if sign > 0 else -result
