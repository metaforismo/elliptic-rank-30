#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define NVARS 8
#define NEQS 6
#define MAXDEG 24

typedef struct {
    int v;
    int d[NVARS];
} Dual;

typedef struct {
    long long jet_solutions;
    long long exact_solutions;
    long long squarefree_solutions;
    long long smooth_rank6_solutions;
    long long smooth_squarefree_solutions;
    long long jacobian_rank_counts[NEQS + 1];
    long long residual_gcd_degree_counts[6];
    int have_example;
    int example_vars[NVARS];
    int example_c4[9];
    int example_c6[13];
    int example_residual[6];
    int example_jacobian_rank;
} SignStats;

static int P;

static inline int modp(int64_t x) {
    x %= P;
    if (x < 0) x += P;
    return (int)x;
}

static inline int addp(int a, int b) {
    int x = a + b;
    if (x >= P) x -= P;
    return x;
}

static inline int subp(int a, int b) {
    int x = a - b;
    if (x < 0) x += P;
    return x;
}

static inline int mulp(int a, int b) {
    return (int)((int64_t)a * b % P);
}

static int powp(int a, int e) {
    int result = 1;
    while (e > 0) {
        if (e & 1) result = mulp(result, a);
        a = mulp(a, a);
        e >>= 1;
    }
    return result;
}

static int invp(int a) {
    if (a == 0) {
        fprintf(stderr, "attempted inversion of zero modulo %d\n", P);
        exit(3);
    }
    return powp(a, P - 2);
}

static int is_prime(int n) {
    if (n < 2) return 0;
    if (n % 2 == 0) return n == 2;
    for (int d = 3; (int64_t)d * d <= n; d += 2) {
        if (n % d == 0) return 0;
    }
    return 1;
}

static int binom_mod(int n, int k) {
    if (k < 0 || k > n) return 0;
    int result = 1;
    for (int i = 1; i <= k; ++i) {
        result = mulp(result, modp(n - k + i));
        result = mulp(result, invp(modp(i)));
    }
    return result;
}

static Dual dconst(int value) {
    Dual result;
    result.v = modp(value);
    memset(result.d, 0, sizeof(result.d));
    return result;
}

static Dual dvar(int value, int index) {
    Dual result = dconst(value);
    result.d[index] = 1;
    return result;
}

static Dual dadd(Dual left, Dual right) {
    Dual result;
    result.v = addp(left.v, right.v);
    for (int i = 0; i < NVARS; ++i) {
        result.d[i] = addp(left.d[i], right.d[i]);
    }
    return result;
}

static Dual dsub(Dual left, Dual right) {
    Dual result;
    result.v = subp(left.v, right.v);
    for (int i = 0; i < NVARS; ++i) {
        result.d[i] = subp(left.d[i], right.d[i]);
    }
    return result;
}

static Dual dscale(Dual value, int scalar) {
    Dual result;
    scalar = modp(scalar);
    result.v = mulp(value.v, scalar);
    for (int i = 0; i < NVARS; ++i) {
        result.d[i] = mulp(value.d[i], scalar);
    }
    return result;
}

static Dual dmul(Dual left, Dual right) {
    Dual result;
    result.v = mulp(left.v, right.v);
    for (int i = 0; i < NVARS; ++i) {
        result.d[i] = addp(
            mulp(left.d[i], right.v),
            mulp(left.v, right.d[i])
        );
    }
    return result;
}

static Dual dsquare(Dual value) {
    return dmul(value, value);
}

static Dual dcube(Dual value) {
    return dmul(dsquare(value), value);
}

static void dacc(Dual *accumulator, int coefficient, Dual term) {
    *accumulator = dadd(*accumulator, dscale(term, coefficient));
}

static void conv_int(
    const int *left,
    int degree_left,
    const int *right,
    int degree_right,
    int *out
) {
    memset(out, 0, (degree_left + degree_right + 1) * sizeof(int));
    for (int i = 0; i <= degree_left; ++i) {
        for (int j = 0; j <= degree_right; ++j) {
            out[i + j] = modp(
                out[i + j] + (int64_t)left[i] * right[j]
            );
        }
    }
}

static void conv_dual(
    const Dual *left,
    int degree_left,
    const Dual *right,
    int degree_right,
    Dual *out
) {
    for (int i = 0; i <= degree_left + degree_right; ++i) {
        out[i] = dconst(0);
    }
    for (int i = 0; i <= degree_left; ++i) {
        for (int j = 0; j <= degree_right; ++j) {
            out[i + j] = dadd(out[i + j], dmul(left[i], right[j]));
        }
    }
}

static void build_c4_dual(const Dual x[NVARS], Dual a[9]) {
    Dual p0 = x[0], p1 = x[1], p2 = x[2], p3 = x[3];
    Dual q0 = x[4], q1 = x[5], q2 = x[6], s = x[7];

    a[0] = dsquare(p0);
    a[1] = dscale(dmul(p0, p1), 2);
    a[2] = dadd(dscale(dmul(p0, p2), 2), dsquare(p1));
    a[3] = dadd(
        dscale(dmul(p0, p3), 2),
        dscale(dmul(p1, p2), 2)
    );

    a[4] = dconst(-3);
    dacc(&a[4], -15, dsquare(p0));
    dacc(&a[4], -20, dmul(p0, p1));
    dacc(&a[4], -12, dmul(p0, p2));
    dacc(&a[4], -6, dmul(p0, p3));
    dacc(&a[4], -6, dsquare(p1));
    dacc(&a[4], -6, dmul(p1, p2));
    dacc(&a[4], 15, dsquare(q0));
    dacc(&a[4], -10, dmul(q0, q1));
    dacc(&a[4], 2, dmul(q0, q2));
    dacc(&a[4], 1, dsquare(q1));
    dacc(&a[4], -1, s);

    a[5] = dconst(8);
    dacc(&a[5], 24, dsquare(p0));
    dacc(&a[5], 30, dmul(p0, p1));
    dacc(&a[5], 16, dmul(p0, p2));
    dacc(&a[5], 6, dmul(p0, p3));
    dacc(&a[5], 8, dsquare(p1));
    dacc(&a[5], 6, dmul(p1, p2));
    dacc(&a[5], -24, dsquare(q0));
    dacc(&a[5], 18, dmul(q0, q1));
    dacc(&a[5], -4, dmul(q0, q2));
    dacc(&a[5], -2, dsquare(q1));
    dacc(&a[5], 3, s);

    a[6] = dconst(-6);
    dacc(&a[6], -10, dsquare(p0));
    dacc(&a[6], -12, dmul(p0, p1));
    dacc(&a[6], -6, dmul(p0, p2));
    dacc(&a[6], -2, dmul(p0, p3));
    dacc(&a[6], -3, dsquare(p1));
    dacc(&a[6], -2, dmul(p1, p2));
    dacc(&a[6], 10, dsquare(q0));
    dacc(&a[6], -8, dmul(q0, q1));
    dacc(&a[6], 2, dmul(q0, q2));
    dacc(&a[6], 1, dsquare(q1));
    dacc(&a[6], -3, s);

    a[7] = s;
    a[8] = dconst(1);
}

static void build_c6_dual(
    const Dual x[NVARS],
    const Dual a[9],
    int e0,
    Dual b[13]
) {
    Dual reversed[13], squared[25], cubed[37], g[13];
    for (int i = 0; i < 13; ++i) reversed[i] = dconst(0);
    for (int k = 0; k <= 8; ++k) reversed[k] = a[8 - k];
    conv_dual(reversed, 8, reversed, 8, squared);
    conv_dual(squared, 16, reversed, 8, cubed);

    int inv2 = invp(2);
    for (int i = 0; i < 13; ++i) g[i] = dconst(0);
    g[0] = dconst(1);
    for (int k = 1; k <= 11; ++k) {
        Dual sum = dconst(0);
        for (int i = 1; i < k; ++i) {
            sum = dadd(sum, dmul(g[i], g[k - i]));
        }
        g[k] = dscale(dsub(cubed[k], sum), inv2);
    }

    for (int i = 0; i < 13; ++i) b[i] = dconst(0);
    b[12] = dconst(1);
    for (int i = 1; i <= 11; ++i) b[i] = g[12 - i];
    b[0] = dscale(dcube(x[0]), e0);
}

static void build_equations(
    const int values[NVARS],
    int e0,
    int e1,
    Dual equations[NEQS],
    int c4[9],
    int c6[13]
) {
    Dual x[NVARS], a[9], b[13];
    for (int i = 0; i < NVARS; ++i) x[i] = dvar(values[i], i);
    build_c4_dual(x, a);
    build_c6_dual(x, a, e0, b);
    for (int i = 0; i <= 8; ++i) c4[i] = a[i].v;
    for (int i = 0; i <= 12; ++i) c6[i] = b[i].v;

    Dual p0 = x[0], p1 = x[1], p2 = x[2], p3 = x[3];
    Dual q0 = x[4], q1 = x[5], q2 = x[6];

    Dual target0 = dscale(dmul(dsquare(p0), p1), 3 * e0);
    Dual target1 = dconst(0);
    dacc(&target1, 3 * e0, dmul(dsquare(p0), p2));
    dacc(&target1, 3 * e0, dmul(p0, dsquare(p1)));
    Dual target2 = dconst(0);
    dacc(&target2, 3 * e0, dmul(dsquare(p0), p3));
    dacc(&target2, 6 * e0, dmul(p0, dmul(p1, p2)));
    dacc(&target2, e0, dcube(p1));
    equations[0] = dsub(b[1], target0);
    equations[1] = dsub(b[2], target1);
    equations[2] = dsub(b[3], target2);

    Dual value0 = dconst(0), value1 = dconst(0), value2 = dconst(0);
    for (int i = 0; i <= 12; ++i) {
        value0 = dadd(value0, b[i]);
        value1 = dadd(value1, dscale(b[i], i));
        value2 = dadd(value2, dscale(b[i], binom_mod(i, 2)));
    }
    Dual target3 = dscale(dcube(q0), e1);
    Dual target4 = dscale(dmul(dsquare(q0), q1), 3 * e1);
    Dual target5 = dconst(0);
    dacc(&target5, 3 * e1, dmul(dsquare(q0), q2));
    dacc(&target5, 3 * e1, dmul(q0, dsquare(q1)));
    equations[3] = dsub(value0, target3);
    equations[4] = dsub(value1, target4);
    equations[5] = dsub(value2, target5);
}

static int matrix_rank_6x8(const Dual equations[NEQS]) {
    int matrix[NEQS][NVARS];
    for (int row = 0; row < NEQS; ++row) {
        for (int column = 0; column < NVARS; ++column) {
            matrix[row][column] = equations[row].d[column];
        }
    }

    int rank = 0;
    for (int column = 0; column < NVARS && rank < NEQS; ++column) {
        int pivot = -1;
        for (int row = rank; row < NEQS; ++row) {
            if (matrix[row][column] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot < 0) continue;
        if (pivot != rank) {
            for (int j = column; j < NVARS; ++j) {
                int temporary = matrix[pivot][j];
                matrix[pivot][j] = matrix[rank][j];
                matrix[rank][j] = temporary;
            }
        }
        int inverse = invp(matrix[rank][column]);
        for (int j = column; j < NVARS; ++j) {
            matrix[rank][j] = mulp(matrix[rank][j], inverse);
        }
        for (int row = 0; row < NEQS; ++row) {
            if (row == rank || matrix[row][column] == 0) continue;
            int factor = matrix[row][column];
            for (int j = column; j < NVARS; ++j) {
                matrix[row][j] = subp(
                    matrix[row][j],
                    mulp(factor, matrix[rank][j])
                );
            }
        }
        ++rank;
    }
    return rank;
}

static int polynomial_degree(const int *polynomial, int maximum_degree) {
    while (maximum_degree >= 0 && polynomial[maximum_degree] == 0) {
        --maximum_degree;
    }
    return maximum_degree;
}

static int polynomial_gcd_degree(
    const int *left_input,
    int left_degree,
    const int *right_input,
    int right_degree
) {
    int left[MAXDEG + 1] = {0};
    int right[MAXDEG + 1] = {0};
    int remainder[MAXDEG + 1] = {0};
    memcpy(left, left_input, (left_degree + 1) * sizeof(int));
    memcpy(right, right_input, (right_degree + 1) * sizeof(int));
    left_degree = polynomial_degree(left, left_degree);
    right_degree = polynomial_degree(right, right_degree);

    while (right_degree >= 0) {
        memcpy(remainder, left, sizeof(remainder));
        int remainder_degree = left_degree;
        int inverse_leading = invp(right[right_degree]);
        while (remainder_degree >= right_degree) {
            int shift = remainder_degree - right_degree;
            int quotient = mulp(
                remainder[remainder_degree], inverse_leading
            );
            for (int i = 0; i <= right_degree; ++i) {
                remainder[i + shift] = subp(
                    remainder[i + shift],
                    mulp(quotient, right[i])
                );
            }
            remainder_degree = polynomial_degree(
                remainder, remainder_degree
            );
        }
        memcpy(left, right, sizeof(left));
        left_degree = right_degree;
        memcpy(right, remainder, sizeof(right));
        right_degree = remainder_degree;
    }
    return left_degree;
}

static int residual_quintic_and_gcd_degree(
    const int discriminant[25],
    int residual[6]
) {
    int quotient_t4[9] = {0};
    for (int i = 0; i <= 8; ++i) quotient_t4[i] = discriminant[i + 4];

    int remainder[9] = {0};
    memcpy(remainder, quotient_t4, sizeof(remainder));
    const int divisor[4] = {
        modp(-1), modp(3), modp(-3), 1
    };
    memset(residual, 0, 6 * sizeof(int));
    for (int degree = 8; degree >= 3; --degree) {
        int coefficient = remainder[degree];
        residual[degree - 3] = coefficient;
        for (int i = 0; i <= 3; ++i) {
            remainder[i + degree - 3] = subp(
                remainder[i + degree - 3],
                mulp(coefficient, divisor[i])
            );
        }
    }
    if (remainder[0] || remainder[1] || remainder[2]) {
        fprintf(stderr, "internal exact division failure modulo %d\n", P);
        exit(4);
    }
    if (residual[5] == 0) {
        fprintf(stderr, "residual polynomial unexpectedly has degree <5\n");
        exit(4);
    }

    int derivative[5] = {0};
    for (int i = 1; i <= 5; ++i) {
        derivative[i - 1] = mulp(i, residual[i]);
    }
    return polynomial_gcd_degree(residual, 5, derivative, 4);
}

static int exact_fiber_orders(
    const int c4[9],
    const int c6[13],
    int discriminant[25]
) {
    int c4_squared[17] = {0};
    int c4_cubed[25] = {0};
    int c6_squared[25] = {0};
    conv_int(c4, 8, c4, 8, c4_squared);
    conv_int(c4_squared, 16, c4, 8, c4_cubed);
    conv_int(c6, 12, c6, 12, c6_squared);
    for (int i = 0; i <= 24; ++i) {
        discriminant[i] = subp(c4_cubed[i], c6_squared[i]);
    }

    for (int i = 0; i < 4; ++i) {
        if (discriminant[i] != 0) return 0;
    }
    if (discriminant[4] == 0) return 0;

    for (int derivative_order = 0; derivative_order < 3; ++derivative_order) {
        int value = 0;
        for (int i = derivative_order; i <= 24; ++i) {
            value = modp(
                value
                + (int64_t)binom_mod(i, derivative_order)
                    * discriminant[i]
            );
        }
        if (value != 0) return 0;
    }
    int third_value = 0;
    for (int i = 3; i <= 24; ++i) {
        third_value = modp(
            third_value
            + (int64_t)binom_mod(i, 3) * discriminant[i]
        );
    }
    if (third_value == 0) return 0;

    for (int i = 13; i <= 24; ++i) {
        if (discriminant[i] != 0) return 0;
    }
    if (discriminant[12] == 0) return 0;
    return 1;
}


static void build_c4_int(const int values[NVARS], int a[9]) {
    int p0 = values[0], p1 = values[1], p2 = values[2], p3 = values[3];
    int q0 = values[4], q1 = values[5], q2 = values[6], s = values[7];
    a[0] = mulp(p0, p0);
    a[1] = modp(2LL * p0 * p1);
    a[2] = modp(2LL * p0 * p2 + (int64_t)p1 * p1);
    a[3] = modp(2LL * p0 * p3 + 2LL * p1 * p2);
    a[4] = modp(
        -15LL * p0 * p0 - 20LL * p0 * p1 - 12LL * p0 * p2
        - 6LL * p0 * p3 - 6LL * p1 * p1 - 6LL * p1 * p2
        + 15LL * q0 * q0 - 10LL * q0 * q1 + 2LL * q0 * q2
        + (int64_t)q1 * q1 - s - 3
    );
    a[5] = modp(
        24LL * p0 * p0 + 30LL * p0 * p1 + 16LL * p0 * p2
        + 6LL * p0 * p3 + 8LL * p1 * p1 + 6LL * p1 * p2
        - 24LL * q0 * q0 + 18LL * q0 * q1 - 4LL * q0 * q2
        - 2LL * q1 * q1 + 3LL * s + 8
    );
    a[6] = modp(
        -10LL * p0 * p0 - 12LL * p0 * p1 - 6LL * p0 * p2
        - 2LL * p0 * p3 - 3LL * p1 * p1 - 2LL * p1 * p2
        + 10LL * q0 * q0 - 8LL * q0 * q1 + 2LL * q0 * q2
        + (int64_t)q1 * q1 - 3LL * s - 6
    );
    a[7] = s;
    a[8] = 1;
}

static void build_c6_base_int(const int c4[9], int c6_base[13]) {
    int reversed[13] = {0};
    int squared[25] = {0};
    int cubed[37] = {0};
    int square_root[13] = {0};
    for (int k = 0; k <= 8; ++k) reversed[k] = c4[8 - k];
    conv_int(reversed, 8, reversed, 8, squared);
    conv_int(squared, 16, reversed, 8, cubed);
    square_root[0] = 1;
    int inverse_two = invp(2);
    for (int k = 1; k <= 11; ++k) {
        int convolution_sum = 0;
        for (int i = 1; i < k; ++i) {
            convolution_sum = modp(
                convolution_sum
                + (int64_t)square_root[i] * square_root[k - i]
            );
        }
        square_root[k] = mulp(
            subp(cubed[k], convolution_sum), inverse_two
        );
    }
    memset(c6_base, 0, 13 * sizeof(int));
    c6_base[12] = 1;
    for (int i = 1; i <= 11; ++i) {
        c6_base[i] = square_root[12 - i];
    }
}

static void print_int_array(const int *values, int length) {
    putchar('[');
    for (int i = 0; i < length; ++i) {
        if (i) printf(", ");
        printf("%d", values[i]);
    }
    putchar(']');
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s prime\n", argv[0]);
        return 2;
    }
    P = atoi(argv[1]);
    if (P < 5 || !is_prime(P)) {
        fprintf(stderr, "prime must be a prime integer >= 5\n");
        return 2;
    }

    SignStats statistics[2][2];
    memset(statistics, 0, sizeof(statistics));
    long long raw_parameter_tuples = 0;
    long long low_filter_sign_attempts = 0;

    for (int p0 = 1; p0 < P; ++p0)
    for (int p1 = 0; p1 < P; ++p1)
    for (int p2 = 0; p2 < P; ++p2)
    for (int p3 = 0; p3 < P; ++p3)
    for (int q0 = 1; q0 < P; ++q0)
    for (int q1 = 0; q1 < P; ++q1)
    for (int q2 = 0; q2 < P; ++q2)
    for (int s = 0; s < P; ++s) {
        ++raw_parameter_tuples;
        int values[NVARS] = {p0, p1, p2, p3, q0, q1, q2, s};
        int c4[9], c6_base[13];
        build_c4_int(values, c4);
        build_c6_base_int(c4, c6_base);

        for (int sign0 = 0; sign0 < 2; ++sign0) {
            int e0 = sign0 ? P - 1 : 1;
            int target_low1 = mulp(
                e0, modp(3LL * p0 * p0 % P * p1)
            );
            int target_low2 = mulp(
                e0,
                modp(
                    3LL * p0 * p0 % P * p2
                    + 3LL * p0 * p1 % P * p1
                )
            );
            int target_low3 = mulp(
                e0,
                modp(
                    3LL * p0 * p0 % P * p3
                    + 6LL * p0 * p1 % P * p2
                    + (int64_t)p1 * p1 % P * p1
                )
            );
            if (
                c6_base[1] != target_low1
                || c6_base[2] != target_low2
                || c6_base[3] != target_low3
            ) {
                continue;
            }

            int c6[13];
            memcpy(c6, c6_base, sizeof(c6));
            c6[0] = mulp(e0, mulp(p0, mulp(p0, p0)));
            int value_at_one = 0;
            int first_jet_at_one = 0;
            int second_jet_at_one = 0;
            for (int i = 0; i <= 12; ++i) {
                value_at_one = addp(value_at_one, c6[i]);
                first_jet_at_one = modp(
                    first_jet_at_one + (int64_t)(i % P) * c6[i]
                );
                second_jet_at_one = modp(
                    second_jet_at_one
                    + (int64_t)binom_mod(i, 2) * c6[i]
                );
            }

            for (int sign1 = 0; sign1 < 2; ++sign1) {
                int e1 = sign1 ? P - 1 : 1;
                ++low_filter_sign_attempts;
                int target_one0 = mulp(e1, mulp(q0, mulp(q0, q0)));
                int target_one1 = mulp(
                    e1, modp(3LL * q0 * q0 % P * q1)
                );
                int target_one2 = mulp(
                    e1,
                    modp(
                        3LL * q0 * q0 % P * q2
                        + 3LL * q0 * q1 % P * q1
                    )
                );
                if (
                    value_at_one != target_one0
                    || first_jet_at_one != target_one1
                    || second_jet_at_one != target_one2
                ) {
                    continue;
                }

                SignStats *entry = &statistics[sign0][sign1];
                ++entry->jet_solutions;
                int discriminant[25] = {0};
                if (!exact_fiber_orders(c4, c6, discriminant)) continue;
                ++entry->exact_solutions;

                Dual equations[NEQS];
                int derivative_c4[9], derivative_c6[13];
                build_equations(
                    values,
                    e0,
                    e1,
                    equations,
                    derivative_c4,
                    derivative_c6
                );
                for (int equation = 0; equation < NEQS; ++equation) {
                    if (equations[equation].v != 0) {
                        fprintf(stderr, "automatic differentiation value mismatch\n");
                        return 4;
                    }
                }
                if (
                    memcmp(c4, derivative_c4, sizeof(c4)) != 0
                    || memcmp(c6, derivative_c6, sizeof(c6)) != 0
                ) {
                    fprintf(stderr, "dual and integer curve construction disagree\n");
                    return 4;
                }

                int jacobian_rank = matrix_rank_6x8(equations);
                ++entry->jacobian_rank_counts[jacobian_rank];
                if (jacobian_rank == 6) ++entry->smooth_rank6_solutions;

                int residual[6] = {0};
                int gcd_degree = residual_quintic_and_gcd_degree(
                    discriminant, residual
                );
                if (gcd_degree < 0 || gcd_degree > 5) {
                    fprintf(stderr, "invalid residual gcd degree\n");
                    return 4;
                }
                ++entry->residual_gcd_degree_counts[gcd_degree];
                if (gcd_degree == 0) {
                    ++entry->squarefree_solutions;
                    if (jacobian_rank == 6) {
                        ++entry->smooth_squarefree_solutions;
                    }
                    if (!entry->have_example) {
                        entry->have_example = 1;
                        memcpy(entry->example_vars, values, sizeof(values));
                        memcpy(entry->example_c4, c4, sizeof(c4));
                        memcpy(entry->example_c6, c6, sizeof(c6));
                        memcpy(
                            entry->example_residual,
                            residual,
                            sizeof(residual)
                        );
                        entry->example_jacobian_rank = jacobian_rank;
                    }
                }
            }
        }
    }

    printf("{\n");
    printf("  \"schema_version\": 1,\n");
    printf("  \"prime\": %d,\n", P);
    printf("  \"parameter_variables\": [\"p0\", \"p1\", \"p2\", \"p3\", \"q0\", \"q1\", \"q2\", \"s\"],\n");
    printf("  \"nonvanishing_constraints\": [\"p0 != 0\", \"q0 != 0\"],\n");
    printf("  \"raw_parameter_tuples\": %lld,\n", raw_parameter_tuples);
    printf("  \"low_filter_sign_attempts\": %lld,\n", low_filter_sign_attempts);
    printf("  \"fiber_orders\": {\"t=0\": 4, \"t=1\": 3, \"t=infinity\": 12},\n");
    printf("  \"expected_dimension_at_rank6_points\": 2,\n");
    printf("  \"sign_data\": [\n");
    for (int sign0 = 0; sign0 < 2; ++sign0) {
        for (int sign1 = 0; sign1 < 2; ++sign1) {
            SignStats *entry = &statistics[sign0][sign1];
            printf("    {\n");
            printf("      \"e0\": %d,\n", sign0 ? -1 : 1);
            printf("      \"e1\": %d,\n", sign1 ? -1 : 1);
            printf("      \"jet_solutions\": %lld,\n", entry->jet_solutions);
            printf("      \"exact_fixed_order_solutions\": %lld,\n", entry->exact_solutions);
            printf("      \"formal_jacobian_rank_counts\": [");
            for (int rank = 0; rank <= NEQS; ++rank) {
                if (rank) printf(", ");
                printf("%lld", entry->jacobian_rank_counts[rank]);
            }
            printf("],\n");
            printf("      \"smooth_rank6_solutions\": %lld,\n", entry->smooth_rank6_solutions);
            printf("      \"residual_gcd_degree_counts\": [");
            for (int degree = 0; degree <= 5; ++degree) {
                if (degree) printf(", ");
                printf("%lld", entry->residual_gcd_degree_counts[degree]);
            }
            printf("],\n");
            printf("      \"squarefree_semistable_solutions\": %lld,\n", entry->squarefree_solutions);
            printf("      \"smooth_squarefree_semistable_solutions\": %lld,\n", entry->smooth_squarefree_solutions);
            printf("      \"example_smooth_semistable\": ");
            if (!entry->have_example) {
                printf("null\n");
            } else {
                printf("{\n");
                printf("        \"parameters\": ");
                print_int_array(entry->example_vars, NVARS);
                printf(",\n");
                printf("        \"c4_coefficients_ascending\": ");
                print_int_array(entry->example_c4, 9);
                printf(",\n");
                printf("        \"c6_coefficients_ascending\": ");
                print_int_array(entry->example_c6, 13);
                printf(",\n");
                printf("        \"residual_quintic_coefficients_ascending\": ");
                print_int_array(entry->example_residual, 6);
                printf(",\n");
                printf("        \"formal_jacobian_rank\": %d\n", entry->example_jacobian_rank);
                printf("      }\n");
            }
            printf("    }%s\n", (sign0 == 1 && sign1 == 1) ? "" : ",");
        }
    }
    printf("  ]\n");
    printf("}\n");
    return 0;
}
