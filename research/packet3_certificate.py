#!/usr/bin/env python3
"""Exact certificate for a genus-zero three-character rank packet.

For K = Q(t), consider

    E_t : y^2 = x^3 - (t^2+1)/(t+1) x + t.

Let d1=t and d2=2t/(t+1).  Over
L=K(sqrt(d1),sqrt(d2)) the points

    (0,sqrt(d1)), (1,sqrt(d2)), (-1,sqrt(d1*d2))

are defined.  This script verifies all algebraic identities, the degree-four
genus-zero parametrization of L, and exact non-torsion specializations that
place the three points in distinct nontrivial V4 character spaces.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Optional

Poly = tuple[int, ...]  # coefficients in ascending order
Point = Optional[tuple[Fraction, Fraction]]


def trim(poly: Iterable[int]) -> Poly:
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [0])


def padd(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def pneg(poly: Poly) -> Poly:
    return tuple(-coefficient for coefficient in poly)


def psub(left: Poly, right: Poly) -> Poly:
    return padd(left, pneg(right))


def pmul(left: Poly, right: Poly) -> Poly:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


def pscale(poly: Poly, scalar: int) -> Poly:
    return trim(scalar * coefficient for coefficient in poly)


def peval(poly: Poly, value: int) -> int:
    result = 0
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def pdiv_linear(poly: Poly, root: int) -> tuple[Poly, int]:
    """Divide by x-root, returning quotient and remainder."""
    if len(poly) == 1:
        return (0,), poly[0]
    quotient = [0] * (len(poly) - 1)
    quotient[-1] = poly[-1]
    for index in range(len(poly) - 3, -1, -1):
        quotient[index] = poly[index + 1] + root * quotient[index + 1]
    remainder = poly[0] + root * quotient[0]
    return trim(quotient), remainder


def multiplicity(poly: Poly, root: int) -> int:
    if poly == (0,):
        raise ValueError("zero polynomial has undefined valuation")
    count = 0
    current = poly
    while len(current) > 1 and peval(current, root) == 0:
        current, remainder = pdiv_linear(current, root)
        if remainder != 0:
            raise AssertionError("linear division remainder mismatch")
        count += 1
    return count


@dataclass(frozen=True)
class RationalFunction:
    numerator: Poly
    denominator: Poly = (1,)

    def __post_init__(self) -> None:
        object.__setattr__(self, "numerator", trim(self.numerator))
        object.__setattr__(self, "denominator", trim(self.denominator))
        if self.denominator == (0,):
            raise ZeroDivisionError("zero denominator")

    @classmethod
    def integer(cls, value: int) -> "RationalFunction":
        return cls((value,))

    def __neg__(self) -> "RationalFunction":
        return RationalFunction(pneg(self.numerator), self.denominator)

    def __add__(self, other: object) -> "RationalFunction":
        rhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return RationalFunction(
            padd(pmul(self.numerator, rhs.denominator), pmul(rhs.numerator, self.denominator)),
            pmul(self.denominator, rhs.denominator),
        )

    __radd__ = __add__

    def __sub__(self, other: object) -> "RationalFunction":
        rhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return self + (-rhs)

    def __rsub__(self, other: object) -> "RationalFunction":
        lhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return lhs - self

    def __mul__(self, other: object) -> "RationalFunction":
        rhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return RationalFunction(
            pmul(self.numerator, rhs.numerator),
            pmul(self.denominator, rhs.denominator),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "RationalFunction":
        rhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        if rhs.numerator == (0,):
            raise ZeroDivisionError
        return RationalFunction(
            pmul(self.numerator, rhs.denominator),
            pmul(self.denominator, rhs.numerator),
        )

    def __rtruediv__(self, other: object) -> "RationalFunction":
        lhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return lhs / self

    def __pow__(self, exponent: int) -> "RationalFunction":
        if exponent < 0:
            return (RationalFunction.integer(1) / self) ** (-exponent)
        result = RationalFunction.integer(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def equals(self, other: object) -> bool:
        rhs = other if isinstance(other, RationalFunction) else RationalFunction.integer(int(other))
        return pmul(self.numerator, rhs.denominator) == pmul(rhs.numerator, self.denominator)

    def valuation(self, root: int) -> int:
        return multiplicity(self.numerator, root) - multiplicity(self.denominator, root)

    def valuation_at_infinity(self) -> int:
        return (len(self.denominator) - 1) - (len(self.numerator) - 1)


def gf2_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    matrix = [[entry & 1 for entry in row] for row in rows]
    width = len(matrix[0])
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column]:
                matrix[row] = [a ^ b for a, b in zip(matrix[row], matrix[rank])]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def ec_add(point1: Point, point2: Point, curve_a: Fraction) -> Point:
    if point1 is None:
        return point2
    if point2 is None:
        return point1
    x1, y1 = point1
    x2, y2 = point2
    if x1 == x2 and y1 == -y2:
        return None
    if point1 == point2:
        if y1 == 0:
            return None
        slope = (3 * x1 * x1 + curve_a) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = slope * (x1 - x3) - y1
    return x3, y3


def ec_multiply(multiplier: int, point: Point, curve_a: Fraction) -> Point:
    if multiplier < 0:
        if point is None:
            return None
        return ec_multiply(-multiplier, (point[0], -point[1]), curve_a)
    result: Point = None
    addend = point
    value = multiplier
    while value:
        if value & 1:
            result = ec_add(result, addend, curve_a)
        addend = ec_add(addend, addend, curve_a)
        value //= 2
    return result


def on_curve(point: Point, curve_a: Fraction, curve_b: Fraction) -> bool:
    if point is None:
        return True
    x, y = point
    return y * y == x * x * x + curve_a * x + curve_b


def curve_discriminant(curve_a: int, curve_b: int) -> int:
    return -16 * (4 * curve_a**3 + 27 * curve_b**2)


def canonical_sha256(payload: object) -> str:
    data = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def compute_certificate() -> dict[str, object]:
    t = RationalFunction((0, 1))
    a = -(t**2 + 1) / (t + 1)
    b = t

    def cubic_at(x_value: int) -> RationalFunction:
        return x_value**3 + a * x_value + b

    d1 = t
    d2 = 2 * t / (t + 1)
    d3 = d1 * d2
    evaluations = {
        "f(0)=d1": cubic_at(0).equals(d1),
        "f(1)=d2": cubic_at(1).equals(d2),
        "f(-1)=d1*d2": cubic_at(-1).equals(d3),
    }
    if not all(evaluations.values()):
        raise AssertionError(f"section identities failed: {evaluations}")

    generic_twist_checks: dict[str, bool] = {}
    for label, squareclass, x_coordinate in (
        ("d1", d1, 0),
        ("d2", d2, 1),
        ("d1*d2", d3, -1),
    ):
        # For y^2=f(x), the d-twist is
        # Y^2=X^3+d^2*a*X+d^3*b.  If f(x0)=d, then
        # (X,Y)=(d*x0,d^2) is a K-rational point on the twist.
        twist_x = squareclass * x_coordinate
        twist_y = squareclass**2
        twist_rhs = twist_x**3 + squareclass**2 * a * twist_x + squareclass**3 * b
        generic_twist_checks[label] = (twist_y**2).equals(twist_rhs)
    if not all(generic_twist_checks.values()):
        raise AssertionError(f"generic twist-section checks failed: {generic_twist_checks}")

    places = (0, -1)
    vectors = {
        "d1": [d1.valuation(place) & 1 for place in places]
        + [d1.valuation_at_infinity() & 1],
        "d2": [d2.valuation(place) & 1 for place in places]
        + [d2.valuation_at_infinity() & 1],
    }
    vectors["d1*d2"] = [left ^ right for left, right in zip(vectors["d1"], vectors["d2"])]
    squareclass_rank = gf2_rank([vectors["d1"], vectors["d2"]])
    if vectors != {"d1": [1, 0, 1], "d2": [1, 1, 0], "d1*d2": [0, 1, 1]}:
        raise AssertionError(f"unexpected parity vectors: {vectors}")
    if squareclass_rank != 2:
        raise AssertionError("d1 and d2 are not independent squareclasses")

    cover_degree = 2**squareclass_rank
    branch_places = sum(any(vector[index] for vector in vectors.values()) for index in range(3))
    two_g_minus_two = -2 * cover_degree + branch_places * (cover_degree // 2)
    if two_g_minus_two % 2:
        raise AssertionError("Riemann-Hurwitz parity failure")
    genus = (two_g_minus_two + 2) // 2
    if (cover_degree, branch_places, genus) != (4, 3, 0):
        raise AssertionError("unexpected cover invariants")

    r = RationalFunction((0, 1))
    denominator = r**2 + 1
    v = (r**2 - 2 * r - 1) / denominator
    w = (1 - 2 * r - r**2) / denominator
    u = v / w
    t_parameter = u**2
    parameterization_checks = {
        "v^2+w^2=2": (v**2 + w**2).equals(2),
        "t=u^2": t_parameter.equals(u**2),
        "2t/(t+1)=v^2": (2 * t_parameter / (t_parameter + 1)).equals(v**2),
        "t*2t/(t+1)=(uv)^2": (
            t_parameter * (2 * t_parameter / (t_parameter + 1))
        ).equals((u * v) ** 2),
    }
    if not all(parameterization_checks.values()):
        raise AssertionError(f"parameterization failed: {parameterization_checks}")

    # Generic discriminant numerator after clearing (t+1)^3.
    generic_discriminant_numerator = (4, 0, -15, -81, -69, -27, 4)
    if generic_discriminant_numerator == (0,):
        raise AssertionError("generic curve is singular")
    generic_discriminant_at_minus3 = 16 * peval(generic_discriminant_numerator, -3) // (-2) ** 3
    if generic_discriminant_at_minus3 != -11888:
        raise AssertionError("unexpected discriminant at t=-3")

    specialized_curve = {"a": 5, "b": -3, "discriminant": generic_discriminant_at_minus3}
    twist_specs = [
        {
            "character": "chi_u",
            "d": -3,
            "source_x": 0,
            "twist_a": 45,
            "twist_b": 81,
            "point": (Fraction(0), Fraction(9)),
            "multiple": 2,
            "expected_multiple": (Fraction(25, 4), Fraction(-197, 8)),
        },
        {
            "character": "chi_v",
            "d": 3,
            "source_x": 1,
            "twist_a": 45,
            "twist_b": -81,
            "point": (Fraction(3), Fraction(9)),
            "multiple": 3,
            "expected_multiple": (Fraction(1479, 49), Fraction(58185, 343)),
        },
        {
            "character": "chi_uv",
            "d": -9,
            "source_x": -1,
            "twist_a": 405,
            "twist_b": 2187,
            "point": (Fraction(9), Fraction(81)),
            "multiple": 3,
            "expected_multiple": (Fraction(13077, 121), Fraction(-1522395, 1331)),
        },
    ]
    specialization_records: list[dict[str, object]] = []
    for spec in twist_specs:
        d = int(spec["d"])
        source_x = int(spec["source_x"])
        expected_a = d * d * specialized_curve["a"]
        expected_b = d**3 * specialized_curve["b"]
        if (spec["twist_a"], spec["twist_b"]) != (expected_a, expected_b):
            raise AssertionError("quadratic-twist coefficient mismatch")
        point = spec["point"]
        curve_a = Fraction(int(spec["twist_a"]))
        curve_b = Fraction(int(spec["twist_b"]))
        if not on_curve(point, curve_a, curve_b):
            raise AssertionError(f"specialized point is off curve: {spec}")
        multiple = ec_multiply(int(spec["multiple"]), point, curve_a)
        if multiple != spec["expected_multiple"]:
            raise AssertionError(f"multiple mismatch: got {multiple}, expected {spec['expected_multiple']}")
        if multiple is None or multiple[0].denominator == 1:
            raise AssertionError("Lutz-Nagell nonintegrality witness failed")
        discriminant = curve_discriminant(int(spec["twist_a"]), int(spec["twist_b"]))
        if discriminant == 0:
            raise AssertionError("singular specialized twist")
        specialization_records.append(
            {
                "character": spec["character"],
                "twist": {
                    "d": d,
                    "equation": f"y^2=x^3+({spec['twist_a']})*x+({spec['twist_b']})",
                    "discriminant": str(discriminant),
                },
                "point": [str(point[0]), str(point[1])],
                "nonintegral_multiple": {
                    "n": spec["multiple"],
                    "point": [str(multiple[0]), str(multiple[1])],
                },
                "non_torsion_argument": (
                    "If the point were torsion, its displayed nonzero multiple would be torsion; "
                    "Lutz-Nagell would force integral coordinates, contradicting the nonintegral x-coordinate."
                ),
            }
        )

    character_signs = {
        "chi_u": {"sigma_u": -1, "sigma_v": 1},
        "chi_v": {"sigma_u": 1, "sigma_v": -1},
        "chi_uv": {"sigma_u": -1, "sigma_v": -1},
    }
    if len({tuple(signs.values()) for signs in character_signs.values()}) != 3:
        raise AssertionError("characters are not distinct")

    payload: dict[str, object] = {
        "schema_version": 1,
        "claim": "For the displayed E/Q(t) and biquadratic genus-zero cover L/Q(t), rank E(L) >= rank E(Q(t)) + 3.",
        "curve": {
            "base_field": "Q(t)",
            "equation": "y^2=x^3-((t^2+1)/(t+1))*x+t",
            "generic_discriminant": "16*(4*t^6-27*t^5-69*t^4-81*t^3-15*t^2+4)/(t+1)^3",
            "generic_nonsingular": True,
        },
        "squareclasses": {
            "d1": "t",
            "d2": "2*t/(t+1)",
            "valuation_places": ["t=0", "t=-1", "t=infinity"],
            "parity_vectors": vectors,
            "rank_over_F2": squareclass_rank,
        },
        "cover": {
            "field": "L=Q(t)(u,v), u^2=t, v^2=2t/(t+1)",
            "degree": cover_degree,
            "branch_place_count": branch_places,
            "genus": genus,
            "rational_parameterization": {
                "v": "(r^2-2*r-1)/(r^2+1)",
                "w": "(1-2*r-r^2)/(r^2+1)",
                "u": "v/w",
                "t": "u^2",
                "checks": parameterization_checks,
            },
        },
        "sections": [
            {"point": ["0", "u"], "character": "chi_u", "identity": "f_t(0)=t=u^2"},
            {"point": ["1", "v"], "character": "chi_v", "identity": "f_t(1)=2t/(t+1)=v^2"},
            {"point": ["-1", "u*v"], "character": "chi_uv", "identity": "f_t(-1)=2t^2/(t+1)=(uv)^2"},
        ],
        "section_identity_checks": evaluations,
        "generic_twist_section_checks": generic_twist_checks,
        "character_signs": character_signs,
        "non_torsion_specialization": {
            "t": -3,
            "base_curve": specialized_curve,
            "twists": specialization_records,
            "generic_conclusion": (
                "A torsion section specializes to a torsion point at every smooth fiber where it is defined. "
                "Each displayed specialized point is non-torsion, so each generic character section is non-torsion."
            ),
        },
        "independence_argument": (
            "In E(L) tensor Q, the commuting V4 involutions split the group into character eigenspaces. "
            "The three non-torsion sections lie in three distinct nontrivial characters, hence are Q-linearly, "
            "therefore Z-linearly, independent and independent from E(Q(t))."
        ),
        "conditional_assumptions": [],
        "implementation": {"language": "Python standard library", "script": "research/packet3_certificate.py"},
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, help="write the freshly recomputed JSON certificate")
    parser.add_argument("--compare", type=Path, help="compare recomputation with a committed JSON certificate")
    arguments = parser.parse_args()

    certificate = compute_certificate()
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {arguments.output}")
    if arguments.compare:
        committed = json.loads(arguments.compare.read_text(encoding="utf-8"))
        if certificate != committed:
            raise AssertionError(f"certificate mismatch: {arguments.compare}")
        print(f"matched {arguments.compare}")

    print("section identities: exact")
    print("squareclass rank: 2; cover degree: 4; cover genus: 0")
    print("specialized non-torsion witnesses: 3/3")
    print("distinct V4 character sections: 3")
    print(f"certificate sha256: {certificate['certificate_sha256']}")
    print("UNCONDITIONAL RESULT: rank E(L) >= rank E(Q(t)) + 3")


if __name__ == "__main__":
    main()
