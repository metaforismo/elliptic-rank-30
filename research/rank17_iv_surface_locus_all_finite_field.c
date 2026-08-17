/* Exhaustive geometric surface-locus enumerator before split-square filtering.
 *
 * The auxiliary jet-sign involution has one representative with e0=1, so this
 * program visits exactly (p-1)*p^5 normalized tuples
 *
 *     (p0,p1,p2,p3,r,s),  p0 != 0.
 *
 * It enforces the five I12+I4+IV surface equations and exact discriminant
 * orders 4,4,12, but records all points whether or not the I4 and IV tangent
 * targets are squares in F_p.  The four split categories are reported
 * separately.  No characteristic-zero or section claim is made.
 */

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int p0, p1, p2, p3, r, s;
    int split_i4, split_iv;
} Candidate;

static int mod_i64(int64_t value, int p) {
    int result = (int)(value % p);
    return result < 0 ? result + p : result;
}

static int inverse_mod(int value, int p) {
    int target = mod_i64(value, p);
    for (int candidate = 1; candidate < p; ++candidate) {
        if ((target * candidate) % p == 1) return candidate;
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

static void build_c6(const int c4[9], int p0, int p, int c6[13]) {
    int reversed[9];
    for (int index = 0; index < 9; ++index) reversed[index] = c4[8 - index];

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
    c6[0] = mod_i64((int64_t)p0 * p0 * p0, p);
    c6[12] = 1;
    for (int order = 1; order < 12; ++order) {
        c6[12 - order] = root[order];
    }
}

static int evaluate_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int index = 0; index < length; ++index) {
        result = mod_i64(result + values[index], p);
    }
    return result;
}

static int derivative_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int index = 1; index < length; ++index) {
        result = mod_i64(result + (int64_t)index * values[index], p);
    }
    return result;
}

static int quadratic_coefficient_at_one(const int *values, int length, int p) {
    int result = 0;
    for (int index = 2; index < length; ++index) {
        result = mod_i64(
            result
            + (int64_t)(index * (index - 1) / 2) * values[index],
            p
        );
    }
    return result;
}

static int accept_surface(
    int p0,
    int p1,
    int p2,
    int p3,
    int r,
    int s,
    int p,
    int *split_i4,
    int *split_iv
) {
    int c4[9], c6[13];
    build_c4(p0, p1, p2, p3, r, s, p, c4);
    build_c6(c4, p0, p, c6);

    int target1 = mod_i64(3LL * p0 * p0 * p1, p);
    int target2 = mod_i64(
        3LL * ((int64_t)p0 * p0 * p2 + (int64_t)p0 * p1 * p1), p
    );
    int target3 = mod_i64(
        3LL * p0 * p0 * p3
        + 6LL * p0 * p1 * p2
        + (int64_t)p1 * p1 * p1,
        p
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
    for (int index = 0; index < 25; ++index) {
        delta[index] = mod_i64(c4_cube[index] - c6_square[index], p);
    }
    for (int index = 0; index < 4; ++index) {
        if (delta[index] != 0) {
            fprintf(stderr, "internal I4 jet failure\n");
            exit(3);
        }
    }
    if (delta[4] == 0 || delta[12] == 0) return 0;
    for (int index = 13; index < 25; ++index) {
        if (delta[index] != 0) {
            fprintf(stderr, "internal I12 recurrence failure\n");
            exit(4);
        }
    }

    *split_i4 = has_square_root(-3 * p0, p);
    *split_iv = has_square_root(-2 * b_at_one, p);
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
    uint64_t categories[2][2] = {{0, 0}, {0, 0}};
    for (int p0 = 1; p0 < p; ++p0) {
        for (int p1 = 0; p1 < p; ++p1)
        for (int p2 = 0; p2 < p; ++p2)
        for (int p3 = 0; p3 < p; ++p3)
        for (int r = 0; r < p; ++r)
        for (int s = 0; s < p; ++s) {
            ++visited;
            int split_i4 = 0, split_iv = 0;
            if (accept_surface(
                p0, p1, p2, p3, r, s, p, &split_i4, &split_iv
            )) {
                Candidate value = {
                    p0, p1, p2, p3, r, s, split_i4, split_iv
                };
                append_candidate(&items, &count, &capacity, value);
                ++categories[split_i4][split_iv];
            }
        }
    }

    printf("{\n  \"prime\": %d,\n", p);
    printf("  \"normalization\": \"e0=1 jet-sign representative\",\n");
    printf("  \"visited_parameter_tuples\": %llu,\n", (unsigned long long)visited);
    printf("  \"candidate_count\": %zu,\n", count);
    printf("  \"split_category_counts\": {\n");
    printf("    \"neither\": %llu,\n", (unsigned long long)categories[0][0]);
    printf("    \"iv_only\": %llu,\n", (unsigned long long)categories[0][1]);
    printf("    \"i4_only\": %llu,\n", (unsigned long long)categories[1][0]);
    printf("    \"both\": %llu\n", (unsigned long long)categories[1][1]);
    printf("  },\n  \"tuples\": [\n");
    for (size_t index = 0; index < count; ++index) {
        Candidate value = items[index];
        printf(
            "    [%d, %d, %d, %d, %d, %d, %d, %d]%s\n",
            value.p0, value.p1, value.p2, value.p3, value.r, value.s,
            value.split_i4, value.split_iv,
            index + 1 == count ? "" : ","
        );
    }
    printf("  ]\n}\n");
    free(items);
    return 0;
}
