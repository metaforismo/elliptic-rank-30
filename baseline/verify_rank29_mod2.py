#!/usr/bin/env python3
"""Exact, standard-library certificate for 29 independent rational points.

The curve is the Elkies--Klagsbrun 2024 record example

    y^2 + x*y = x^3 + a4*x + a6.

The proof has two exact parts:

1. direct rational substitution verifies all 29 points;
2. the reductions of the points are F_2-linearly independent in a product
   of local quotients E(F_p)/2E(F_p).

If an integral relation existed, reducing it in every local quotient would
force all coefficients to be even.  The curve has trivial rational torsion
(certified by #E(F_67)=83 and #E(F_71)=75), so division by two may be
repeated; hence every coefficient is divisible by every power of two and
must be zero.

No floating point arithmetic, BSD, GRH, analytic-rank estimate, SageMath,
Magma, or PARI/GP is used by this certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Optional

ROOT = Path(__file__).resolve().parents[1]
CURVE_PATH = ROOT / "curve.json"
POINTS_PATH = ROOT / "points.json"
DEFAULT_CERTIFICATE_PATH = ROOT / "baseline" / "rank29_mod2_certificate.json"

LOCAL_PRIMES = (19, 23, 37, 59, 107, 127, 173, 179, 263, 317, 379, 389, 397, 419, 463)
TORSION_PRIMES = (67, 71)
EXPECTED_LOCAL_ORDERS = {
    19: 28,
    23: 32,
    37: 48,
    59: 72,
    107: 128,
    127: 136,
    173: 200,
    179: 200,
    263: 288,
    317: 348,
    379: 400,
    389: 416,
    397: 416,
    419: 432,
    463: 496,
}
EXPECTED_TORSION_ORDERS = {67: 83, 71: 75}
EXPECTED_MATRIX_RANK = 29

Point = Optional[tuple[int, int]]
RationalPoint = tuple[Fraction, Fraction]


def load_inputs() -> tuple[int, int, int, int, int, list[RationalPoint]]:
    curve = json.loads(CURVE_PATH.read_text(encoding="utf-8"))
    ainvs = [int(value) for value in curve["model"]["a_invariants"]]
    if ainvs[:3] != [1, 0, 0]:
        raise AssertionError(f"unexpected a1,a2,a3: {ainvs[:3]}")
    a1, a2, a3, a4, a6 = ainvs

    point_data = json.loads(POINTS_PATH.read_text(encoding="utf-8"))
    points: list[RationalPoint] = []
    for expected_index, item in enumerate(point_data["points"], start=1):
        point_id, x_num, x_den, y_num, y_den = item
        if point_id != f"P{expected_index}":
            raise AssertionError(f"unexpected point id: {point_id!r}")
        x = Fraction(int(x_num), int(x_den))
        y = Fraction(int(y_num), int(y_den))
        points.append((x, y))
    if len(points) != 29:
        raise AssertionError(f"expected 29 points, found {len(points)}")
    return a1, a2, a3, a4, a6, points


def generalized_invariants(
    a1: int, a2: int, a3: int, a4: int, a6: int
) -> dict[str, int]:
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6
    discriminant = (
        -b2 * b2 * b8
        - 8 * b4 * b4 * b4
        - 27 * b6 * b6
        + 9 * b2 * b4 * b6
    )
    if c4**3 - c6**2 != 1728 * discriminant:
        raise AssertionError("invariant identity c4^3-c6^2=1728*Delta failed")
    return {
        "b2": b2,
        "b4": b4,
        "b6": b6,
        "b8": b8,
        "c4": c4,
        "c6": c6,
        "discriminant": discriminant,
    }


def verify_rational_points(a4: int, a6: int, points: Iterable[RationalPoint]) -> None:
    for index, (x, y) in enumerate(points, start=1):
        residual = y * y + x * y - (x * x * x + a4 * x + a6)
        if residual != 0:
            raise AssertionError(f"P{index} is not on the curve: residual={residual}")


def short_model_coefficients(a4: int, a6: int) -> tuple[int, int]:
    # X = 36*x + 3, Y = 108*(2*y+x)
    short_a = 1296 * a4 - 27
    short_b = 46656 * a6 - 3888 * a4 + 54
    return short_a, short_b


def inverse_mod(value: int, modulus: int) -> int:
    return pow(value % modulus, -1, modulus)


def point_key(point: Point) -> tuple[int, int, int]:
    if point is None:
        return (-1, -1, -1)
    return (0, point[0], point[1])


def ec_add(point1: Point, point2: Point, prime: int, a: int) -> Point:
    if point1 is None:
        return point2
    if point2 is None:
        return point1
    x1, y1 = point1
    x2, y2 = point2
    if x1 == x2 and (y1 + y2) % prime == 0:
        return None
    if point1 == point2:
        if y1 % prime == 0:
            return None
        slope = (3 * x1 * x1 + a) * inverse_mod(2 * y1, prime) % prime
    else:
        slope = (y2 - y1) * inverse_mod(x2 - x1, prime) % prime
    x3 = (slope * slope - x1 - x2) % prime
    y3 = (slope * (x1 - x3) - y1) % prime
    return (x3, y3)


def enumerate_curve_points(prime: int, a: int, b: int) -> list[Point]:
    square_roots: dict[int, list[int]] = {}
    for y in range(prime):
        square_roots.setdefault(y * y % prime, []).append(y)
    points: list[Point] = [None]
    for x in range(prime):
        rhs = (x * x * x + a * x + b) % prime
        points.extend((x, y) for y in square_roots.get(rhs, ()))
    return points


def quotient_by_doubling(
    group: list[Point], prime: int, a: int
) -> tuple[dict[Point, tuple[int, ...]], list[Point], list[list[Point]]]:
    doubled = {ec_add(point, point, prime, a) for point in group}
    if None not in doubled:
        raise AssertionError("2E(F_p) does not contain the identity")

    unassigned = set(group)
    cosets: list[list[Point]] = []
    point_to_coset: dict[Point, int] = {}
    while unassigned:
        representative = min(unassigned, key=point_key)
        coset_set = {ec_add(representative, element, prime, a) for element in doubled}
        coset = sorted(coset_set, key=point_key)
        index = len(cosets)
        for point in coset:
            if point in point_to_coset:
                raise AssertionError("cosets overlap")
            point_to_coset[point] = index
        unassigned.difference_update(coset_set)
        cosets.append(coset)

    quotient_size = len(cosets)
    if quotient_size not in (1, 2, 4):
        raise AssertionError(f"unexpected quotient size {quotient_size}")
    zero_index = point_to_coset[None]
    nonzero_indices = [index for index in range(quotient_size) if index != zero_index]

    index_to_vector: dict[int, tuple[int, ...]]
    basis_representatives: list[Point]
    if quotient_size == 1:
        index_to_vector = {zero_index: ()}
        basis_representatives = []
    elif quotient_size == 2:
        index_to_vector = {zero_index: (0,), nonzero_indices[0]: (1,)}
        basis_representatives = [cosets[nonzero_indices[0]][0]]
    else:
        first, second, third = nonzero_indices
        rep1 = cosets[first][0]
        rep2 = cosets[second][0]
        sum_index = point_to_coset[ec_add(rep1, rep2, prime, a)]
        if sum_index == zero_index:
            raise AssertionError("chosen quotient generators are dependent")
        if sum_index not in (first, second, third):
            raise AssertionError("invalid quotient addition")
        remaining = next(
            index for index in nonzero_indices if index not in (first, second)
        )
        if sum_index != remaining:
            raise AssertionError("quotient addition table is inconsistent")
        index_to_vector = {
            zero_index: (0, 0),
            first: (1, 0),
            second: (0, 1),
            remaining: (1, 1),
        }
        basis_representatives = [rep1, rep2]

    point_to_vector = {
        point: index_to_vector[point_to_coset[point]] for point in group
    }
    return point_to_vector, basis_representatives, cosets


def reduce_rational_point(point: RationalPoint, prime: int) -> Point:
    x, y = point
    if x.denominator % prime == 0 or y.denominator % prime == 0:
        raise AssertionError(f"point is nonintegral at p={prime}")
    x_mod = x.numerator * inverse_mod(x.denominator, prime) % prime
    y_mod = y.numerator * inverse_mod(y.denominator, prime) % prime
    return ((36 * x_mod + 3) % prime, (108 * (2 * y_mod + x_mod)) % prime)


def gf2_rref(rows: list[list[int]]) -> tuple[int, list[int], list[list[int]]]:
    if not rows:
        return 0, [], []
    width = len(rows[0])
    matrix = [[entry & 1 for entry in row] for row in rows]
    if any(len(row) != width for row in matrix):
        raise AssertionError("ragged binary matrix")

    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        source = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        for row in range(len(matrix)):
            if row != pivot_row and matrix[row][column]:
                matrix[row] = [
                    left ^ right
                    for left, right in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row, pivot_columns, matrix


def point_to_json(point: Point) -> object:
    if point is None:
        return "O"
    return [point[0], point[1]]


def canonical_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compute_certificate() -> dict[str, object]:
    a1, a2, a3, a4, a6, rational_points = load_inputs()
    invariants = generalized_invariants(a1, a2, a3, a4, a6)
    if invariants["discriminant"] == 0:
        raise AssertionError("singular curve")
    verify_rational_points(a4, a6, rational_points)

    short_a, short_b = short_model_coefficients(a4, a6)
    local_records: list[dict[str, object]] = []
    matrix_rows: list[list[int]] = [[] for _ in rational_points]

    for prime in LOCAL_PRIMES:
        if invariants["discriminant"] % prime == 0:
            raise AssertionError(f"p={prime} is not a good-reduction prime")
        finite_a = short_a % prime
        finite_b = short_b % prime
        group = enumerate_curve_points(prime, finite_a, finite_b)
        if len(group) != EXPECTED_LOCAL_ORDERS[prime]:
            raise AssertionError(
                f"#E(F_{prime})={len(group)}, expected {EXPECTED_LOCAL_ORDERS[prime]}"
            )
        vector_map, basis, cosets = quotient_by_doubling(group, prime, finite_a)
        dimension = {1: 0, 2: 1, 4: 2}[len(cosets)]
        if dimension != 2:
            raise AssertionError(
                f"local quotient at p={prime} has dimension {dimension}"
            )

        vectors: list[str] = []
        for index, rational_point in enumerate(rational_points):
            reduced = reduce_rational_point(rational_point, prime)
            if reduced not in vector_map:
                raise AssertionError(f"P{index + 1} does not reduce to E(F_{prime})")
            vector = vector_map[reduced]
            matrix_rows[index].extend(vector)
            vectors.append("".join(str(bit) for bit in vector))

        local_records.append(
            {
                "prime": prime,
                "group_order": len(group),
                "quotient_dimension": dimension,
                "quotient_basis_representatives": [
                    point_to_json(point) for point in basis
                ],
                "point_vectors": vectors,
            }
        )

    rank, pivot_columns, _rref = gf2_rref(matrix_rows)
    if rank != EXPECTED_MATRIX_RANK:
        raise AssertionError(f"local image matrix has rank {rank}, expected 29")

    torsion_records: list[dict[str, int]] = []
    for prime in TORSION_PRIMES:
        if invariants["discriminant"] % prime == 0:
            raise AssertionError(f"torsion prime p={prime} has bad reduction")
        group_order = len(
            enumerate_curve_points(prime, short_a % prime, short_b % prime)
        )
        if group_order != EXPECTED_TORSION_ORDERS[prime]:
            raise AssertionError(
                f"#E(F_{prime})={group_order}, expected {EXPECTED_TORSION_ORDERS[prime]}"
            )
        torsion_records.append({"prime": prime, "group_order": group_order})

    payload: dict[str, object] = {
        "schema_version": 1,
        "candidate_id": "elkies_klagsbrun_2024_rank29",
        "claim": (
            "The 29 listed rational points are Z-linearly independent; "
            "rank E(Q) >= 29."
        ),
        "proof_type": "exact reduction to a product of E(F_p)/2E(F_p)",
        "curve": {
            "a_invariants": [str(a1), str(a2), str(a3), str(a4), str(a6)],
            "discriminant": str(invariants["discriminant"]),
            "short_model": {
                "equation": "Y^2 = X^3 + A*X + B",
                "A": str(short_a),
                "B": str(short_b),
                "map": {"X": "36*x+3", "Y": "108*(2*y+x)"},
            },
        },
        "point_count": len(rational_points),
        "point_membership": "verified exactly over Q",
        "torsion_certificate": {
            "reductions": torsion_records,
            "conclusion": "E(Q)_tors is trivial",
            "argument": (
                "For a rational torsion point of prime order ell, reduction is "
                "injective at every good prime p != ell. If ell is neither 67 "
                "nor 71, then ell divides gcd(83,75)=1. If ell=67, reduction "
                "at 71 would force 67|75; if ell=71, reduction at 67 would force "
                "71|83. Hence no prime-order rational torsion point exists."
            ),
        },
        "local_quotients": local_records,
        "binary_matrix": {
            "orientation": (
                "29 point rows by 30 concatenated local-quotient coordinates"
            ),
            "rows": ["".join(str(bit) for bit in row) for row in matrix_rows],
            "rank": rank,
            "pivot_columns_zero_based": pivot_columns,
        },
        "independence_argument": (
            "Any integral relation reduces in every quotient E(F_p)/2E(F_p). "
            "The displayed binary matrix has row rank 29, so every coefficient "
            "is even. Halving the coefficients gives a point killed by 2; "
            "trivial rational torsion makes it O. Iterating forces every "
            "coefficient to be divisible by every power of 2, hence all are zero."
        ),
        "conditional_assumptions": [],
        "implementation": {
            "language": "Python standard library",
            "script": "baseline/verify_rank29_mod2.py",
        },
    }
    payload["certificate_sha256"] = canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-certificate",
        type=Path,
        metavar="PATH",
        help="write the freshly computed canonical JSON certificate",
    )
    parser.add_argument(
        "--certificate",
        type=Path,
        default=DEFAULT_CERTIFICATE_PATH,
        help="committed certificate to compare against",
    )
    args = parser.parse_args()

    computed = compute_certificate()
    if args.write_certificate is not None:
        args.write_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.write_certificate.write_text(
            json.dumps(computed, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {args.write_certificate}")
    else:
        committed = json.loads(args.certificate.read_text(encoding="utf-8"))
        if committed != computed:
            raise AssertionError(
                f"committed certificate {args.certificate} does not match recomputation"
            )

    matrix = computed["binary_matrix"]
    print("exact point checks: 29/29")
    print(f"nonsingular discriminant: {computed['curve']['discriminant']}")
    print("torsion: trivial (#E(F_67)=83, #E(F_71)=75)")
    print(
        "local quotient matrix: "
        f"{len(matrix['rows'])}x{len(matrix['rows'][0])}, rank={matrix['rank']}"
    )
    print(f"certificate sha256: {computed['certificate_sha256']}")
    print("UNCONDITIONAL RESULT: rank E(Q) >= 29")


if __name__ == "__main__":
    main()
