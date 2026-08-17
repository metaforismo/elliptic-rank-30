/* Exhaustive finite-field enumerator for the normalized I12+I4+IV surface chart.
 *
 * The output contains only the accepted parameter tuples.  A separate Python
 * verifier rebuilds c4, c6, the discriminant, and every local condition from
 * those tuples.  Arithmetic uses int values reduced modulo the requested odd
 * prime; the intended production run is p=11.
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int e0, p0, p1, p2, p3, r, s;
} Candidate;

static int mod_i64(int64_t value, int p) {
    int result = (int)(value % p);
    return result < 0 ? result + p : result;
}

static int inverse_mod(int value, int p) {
    int a = mod_i64(value, p);
    for (int x = 1; x < p; ++x) {
        if ((a * x) % p == 1) return x;
    }
    fprintf(stderr, "noninvertible value %d modulo %d\n", value, p);
    exit(2);
}

static int has_square_root(int value, int p) {
    int target = mod_i64(value, p);
    for (int x = 0; x < p; ++x) {
        if ((x * x) % p == target) return 1;
    }
    return 0;
}

static void convolution(
    const int *left,
    int left_length,
    const int *right,
    int right_length,
    int *output,
    int output_length,
    int p
) {
    memset(output, 0, (size_t)output_length * sizeof(int));
    for (int i = 0; i < left_length; ++i) {
        if (!left[i]) continue;
        for (int j = 0; j < right_length && i + j < output_length; ++j) {
            if (!right[j]) continue;
            output[i + j] = mod_i64(
                output[i + j] + (int64_t)left[i] * right[j], p
            );
        }
    }
}

static void build_c4(
    int p0,
    int p1,
    int p2,
    int p3,
    int r,
    int s,
    int p,
    int c4[9]
) {
    int a0 = mod_i64((int64_t)p0 * p0, p);
    int a1 = mod_i64(2LL * p0 * p1, p);
    int a2 = mod_i64(2LL * p0 * p2 + (int64_t)p1 * p1, p);
    int a3 = mod_i64(2LL * p0 * p3 + 2LL * p1 * p2, p);
    int l0 = mod_i64((int64_t)a0 + a1 + a2 + a3 + s + 1, p);
    int l1 = mod_i64((int64_t)a1 + 2LL * a2 + 3LL * a3 + 7LL * s + 8, p);
    int a4 = mod_i64((int64_t)l1 - 5LL * l0 + r, p);
    int a5 = mod_i64(4LL * l0 - l1 - 2LL * r, p);
    int values[9] = {a0, a1, a2, a3, a4, a5, r, s, 1};
    memcpy(c4, values, sizeof(values));
}

static void build_c6(
    const int c4[9],
    int c6_constant,
    int p,
    int c6[13]
) {
    int reversed[9];
    for (int i = 0; i < 9; ++i) reversed[i] = c4[8 - i];

    int square[12], cube[12];
    convolution(reversed, 9, reversed, 9, square, 12, p);
    convolution(square, 12, reversed, 9, cube, 12, p);

    int root[12] = {0};
    root[0] = 1;
    int inverse_two = inverse_mod(2, p);
    for (int order = 1; order < 12; ++order) {
        int correction = 0;
        for (int index = 1; index < order; ++index) {
            correction = mod_i64(
                correction + (int64_t)root[index] * root[order - index], p
            );
        }
        root[order] = mod_i64(
            (int64_t)(cube[order] - correction) * inverse_two, p
        );
    }

    memset(c6, 0, 13 * sizeof(int));
    c6[0] = mod_i64(c6_constant, p);
    c6[12] = 1;
    for (int order = 1; order < 12; ++order) {
        c6[12 - order] = root[order];
    }
}

static int evaluate_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int i = 0; i < length; ++i) result = mod_i64(result + values[i], p);
    return result;
}

static int derivative_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int i = 1; i < length; ++i) {
        result = mod_i64(result + (int64_t)i * values[i], p);
    }
    return result;
}

static int quadratic_coefficient_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int i = 2; i < length; ++i) {
        result = mod_i64(
            result + (int64_t)(i * (i - 1) / 2) * values[i], p
        );
    }
    return result;
}

static int accept_surface(
    int e0,
    int p0,
    int p1,
    int p2,
    int p3,
    int r,
    int s,
    int p
) {
    int c4[9], c6[13];
    build_c4(p0, p1, p2, p3, r, s, p, c4);
    build_c6(c4, mod_i64((int64_t)e0 * p0 * p0 * p0, p), p, c6);

    int target1 = mod_i64(3LL * e0 * p0 * p0 * p1, p);
    int target2 = mod_i64(
        3LL * e0 * ((int64_t)p0 * p0 * p2 + (int64_t)p0 * p1 * p1), p
    );
    int target3 = mod_i64(
        (int64_t)e0 * (
            3LL * p0 * p0 * p3
            + 6LL * p0 * p1 * p2
            + (int64_t)p1 * p1 * p1
        ), p
    );
    if (c6[1] != target1 || c6[2] != target2 || c6[3] != target3) return 0;
    if (evaluate_at_one(c6, 13, p) != 0) return 0;
    if (derivative_at_one(c6, 13, p) != 0) return 0;

    int b_at_one = quadratic_coefficient_at_one(c6, 13, p);
    if (b_at_one == 0) return 0;

    int c4_square[17], c4_cube[25], c6_square[25], delta[25];
    convolution(c4, 9, c4, 9, c4_square, 17, p);
    convolution(c4_square, 17, c4, 9, c4_cube, 25, p);
    convolution(c6, 13, c6, 13, c6_square, 25, p);
    for (int i = 0; i < 25; ++i) {
        delta[i] = mod_i64(c4_cube[i] - c6_square[i], p);
    }
    for (int i = 0; i < 4; ++i) {
        if (delta[i] != 0) {
            fprintf(stderr, "internal I4 jet failure\n");
            exit(3);
        }
    }
    if (delta[4] == 0 || delta[12] == 0) return 0;
    for (int i = 13; i < 25; ++i) {
        if (delta[i] != 0) {
            fprintf(stderr, "internal I12 recurrence failure\n");
            exit(4);
        }
    }

    if (!has_square_root(-3 * e0 * p0, p)) return 0;
    if (!has_square_root(-2 * b_at_one, p)) return 0;
    return 1;
}

static void append_candidate(
    Candidate **items,
    size_t *count,
    size_t *capacity,
    Candidate value
) {
    if (*count == *capacity) {
        size_t next_capacity = *capacity ? 2 * *capacity : 64;
        Candidate *next = realloc(*items, next_capacity * sizeof(Candidate));
        if (!next) {
            perror("realloc");
            exit(5);
        }
        *items = next;
        *capacity = next_capacity;
    }
    (*items)[(*count)++] = value;
}

int main(int argc, char **argv) {
    int p = 11;
    if (argc >= 2) {
        char *end = NULL;
        errno = 0;
        long parsed = strtol(argv[1], &end, 10);
        if (errno || !end || *end || parsed <= 3 || parsed > 46337) {
            fprintf(stderr, "usage: %s PRIME_GREATER_THAN_3\n", argv[0]);
            return 2;
        }
        p = (int)parsed;
    }

    Candidate *items = NULL;
    size_t count = 0, capacity = 0;
    uint64_t visited = 0;
    int signs[2] = {1, p - 1};
    for (int sign_index = 0; sign_index < 2; ++sign_index) {
        int e0 = signs[sign_index];
        for (int p0 = 1; p0 < p; ++p0) {
            for (int p1 = 0; p1 < p; ++p1)
            for (int p2 = 0; p2 < p; ++p2)
            for (int p3 = 0; p3 < p; ++p3)
            for (int r = 0; r < p; ++r)
            for (int s = 0; s < p; ++s) {
                ++visited;
                if (accept_surface(e0, p0, p1, p2, p3, r, s, p)) {
                    Candidate value = {e0, p0, p1, p2, p3, r, s};
                    append_candidate(&items, &count, &capacity, value);
                }
            }
        }
    }

    printf("{\n  \"prime\": %d,\n", p);
    printf("  \"visited_parameter_tuples\": %llu,\n", (unsigned long long)visited);
    printf("  \"candidate_count\": %zu,\n", count);
    printf("  \"tuples\": [\n");
    for (size_t index = 0; index < count; ++index) {
        Candidate value = items[index];
        printf(
            "    [%d, %d, %d, %d, %d, %d, %d]%s\n",
            value.e0, value.p0, value.p1, value.p2,
            value.p3, value.r, value.s,
            index + 1 == count ? "" : ","
        );
    }
    printf("  ]\n}\n");
    free(items);
    return 0;
}
