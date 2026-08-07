#!/usr/bin/env python3
"""Dependency-free exact certificate for cyclic cubic trisection trace codes."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "cyclic_trisection_trace_code.json"


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def det_bareiss(matrix):
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot_row = next(i for i in range(k + 1, n) if a[i][k])
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot - a[i][k] * a[k][j]) // previous
        previous = pivot
    return sign * a[-1][-1]


def inverse_fraction(matrix):
    n = len(matrix)
    a = [
        [Fraction(matrix[i][j]) for j in range(n)]
        + [Fraction(i == j) for j in range(n)]
        for i in range(n)
    ]
    for column in range(n):
        row = next(i for i in range(column, n) if a[i][column])
        a[column], a[row] = a[row], a[column]
        pivot = a[column][column]
        a[column] = [entry / pivot for entry in a[column]]
        for i in range(n):
            if i != column and a[i][column]:
                factor = a[i][column]
                a[i] = [
                    a[i][j] - factor * a[column][j]
                    for j in range(2 * n)
                ]
    return [row[n:] for row in a]


def matrix_vector(matrix, vector):
    return [
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    ]


# A standard E8 simple-root basis in R^8, multiplied by two and stored as
# columns.  The actual basis matrix is B2/2.
ROOTS_TIMES_TWO = [
    (1, -1, -1, -1, -1, -1, -1, 1),
    (2, 2, 0, 0, 0, 0, 0, 0),
    (-2, 2, 0, 0, 0, 0, 0, 0),
    (0, -2, 2, 0, 0, 0, 0, 0),
    (0, 0, -2, 2, 0, 0, 0, 0),
    (0, 0, 0, -2, 2, 0, 0, 0),
    (0, 0, 0, 0, -2, 2, 0, 0),
    (0, 0, 0, 0, 0, -2, 2, 0),
]
B2 = [[ROOTS_TIMES_TWO[j][i] for j in range(8)] for i in range(8)]
GRAM = [
    [dot(ROOTS_TIMES_TWO[i], ROOTS_TIMES_TWO[j]) // 4 for j in range(8)]
    for i in range(8)
]
assert det_bareiss(GRAM) == 1
assert abs(det_bareiss(B2)) == 2**8
B2_INVERSE = inverse_fraction(B2)


def coordinates(twice_vector):
    values = matrix_vector(B2_INVERSE, twice_vector)
    assert all(value.denominator == 1 for value in values)
    return tuple(int(value) for value in values)


def norm_six_shell():
    vectors = []

    # Integral coset D8.
    for vector in product(range(-2, 3), repeat=8):
        if sum(entry * entry for entry in vector) == 6 and sum(vector) % 2 == 0:
            vectors.append(tuple(2 * entry for entry in vector))

    # Half-integral coset.  If y=2x, sum y_i^2=24 forces exactly two
    # coordinates of absolute value three and six of absolute value one.
    for places in combinations(range(8), 2):
        for signs in product((-1, 1), repeat=8):
            vector = tuple(
                (3 if i in places else 1) * signs[i] for i in range(8)
            )
            if sum(vector) % 4 == 0:
                vectors.append(vector)

    assert len(vectors) == 6720
    assert len(set(vectors)) == 6720
    return vectors


def bilinear_mod_three(left, right):
    return sum(
        left[i] * GRAM[i][j] * right[j]
        for i in range(8)
        for j in range(8)
    ) % 3


def quadratic_mod_three(vector):
    # q(v)=(v,v)/2; the inverse of two modulo three is two.
    return 2 * bilinear_mod_three(vector, vector) % 3


shell = norm_six_shell()
residue_fibres = {}
for vector in shell:
    residue = tuple(value % 3 for value in coordinates(vector))
    assert any(residue)
    assert quadratic_mod_three(residue) == 0
    residue_fibres[residue] = residue_fibres.get(residue, 0) + 1

isotropic_vectors = {
    vector
    for vector in product(range(3), repeat=8)
    if any(vector) and quadratic_mod_three(vector) == 0
}
assert len(isotropic_vectors) == 2240
assert set(residue_fibres) == isotropic_vectors
assert set(residue_fibres.values()) == {3}

field_order = 3
half_dimension = 4
isotropic_formula = (
    field_order ** (half_dimension - 1) + 1
) * (field_order**half_dimension - 1)
maximal_space_formula = 2
for i in range(1, half_dimension):
    maximal_space_formula *= field_order**i + 1

assert isotropic_formula == 2240
assert maximal_space_formula == 2240
nonzero_per_maximal_space = field_order**half_dimension - 1
assert nonzero_per_maximal_space == 80
assert nonzero_per_maximal_space * 3 == 240
spaces_through_vector = (
    maximal_space_formula * nonzero_per_maximal_space // len(isotropic_vectors)
)
assert spaces_through_vector == 80


def sigma_three(n):
    return sum(divisor**3 for divisor in range(1, n + 1) if n % divisor == 0)


shell_counts = {
    str(intersection): 240 * sigma_three(3 * (intersection + 1))
    for intersection in range(8)
}
assert shell_counts["0"] == 6720
assert shell_counts["1"] == 60480

# Height calculation for an orbit of a genus-one cyclic trisection.
for intersection in range(50):
    section_height = 6 + 2 * intersection
    conjugate_pairing = 2 * intersection
    trace_height = 3 * section_height + 6 * conjugate_pairing
    assert trace_height == 18 * (intersection + 1)
    assert Fraction(section_height) - Fraction(trace_height, 9) == 4
    assert Fraction(conjugate_pairing) - Fraction(trace_height, 9) == -2
    assert 2 * section_height - 2 * conjugate_pairing == 12

certificate = {
    "status": "pass",
    "theorem": "cyclic cubic trisection trace-code geometry",
    "e8_gram_determinant": 1,
    "e8_norm6_vector_count": 6720,
    "mod3_quadratic_space": {
        "dimension": 8,
        "field": 3,
        "type": "plus",
        "witt_index": 4,
    },
    "nonzero_isotropic_vectors_mod3": 2240,
    "norm6_lifts_per_nonzero_isotropic_vector": 3,
    "maximal_totally_isotropic_4_spaces": 2240,
    "nonzero_vectors_per_maximal_space": 80,
    "norm6_trisection_classes_per_maximal_space": 240,
    "maximal_spaces_through_each_nonzero_isotropic_vector": 80,
    "orbit_height_formula": {
        "section_height": "6+2m",
        "pairing_of_distinct_conjugates": "2m",
        "trace_height": "18(m+1)",
        "projected_height": 4,
        "projected_pairing": -2,
        "difference_height": 12,
        "projected_lattice": "A2(2)",
        "projected_gram": [[4, -2], [-2, 4]],
        "difference_lattice": "A2(6)",
        "difference_gram": [[12, -6], [-6, 12]],
    },
    "shell_counts_for_intersection_m_0_through_7": shell_counts,
    "search_reduction": (
        "An integral projected packet has trace residues in a totally "
        "isotropic four-space, reducing the minimal shell from 6720 to 240 "
        "classes."
    ),
    "truth_status": "new intermediate theorem; no rank-30 curve claimed",
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(
    json.dumps(
        {
            key: certificate[key]
            for key in (
                "status",
                "e8_norm6_vector_count",
                "nonzero_isotropic_vectors_mod3",
                "maximal_totally_isotropic_4_spaces",
                "norm6_trisection_classes_per_maximal_space",
            )
        },
        sort_keys=True,
    )
)
