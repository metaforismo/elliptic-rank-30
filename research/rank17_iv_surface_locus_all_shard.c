/* Exact p0-sharded wrapper for rank17_iv_surface_locus_all_finite_field.c.
 *
 * Usage:
 *
 *   shard PRIME P0_START_INCLUSIVE P0_END_EXCLUSIVE
 *
 * The range must lie in 1 <= start < end <= PRIME.  All remaining five
 * coordinates are exhausted.  The included implementation supplies the exact
 * invariant construction and surface predicate; only its main function is
 * renamed.  Concatenating a disjoint partition of [1, PRIME) therefore gives
 * the same census as the unsharded program.
 */

#define main rank17_iv_unsharded_main
#include "rank17_iv_surface_locus_all_finite_field.c"
#undef main

static int parse_integer(const char *text, const char *label) {
    char *end = NULL;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno || !end || *end || parsed < 0 || parsed > 46337) {
        fprintf(stderr, "invalid %s: %s\n", label, text);
        exit(2);
    }
    return (int)parsed;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(
            stderr,
            "usage: %s PRIME P0_START_INCLUSIVE P0_END_EXCLUSIVE\n",
            argv[0]
        );
        return 2;
    }

    int p = parse_integer(argv[1], "prime");
    int start = parse_integer(argv[2], "p0 start");
    int end = parse_integer(argv[3], "p0 end");
    if (p <= 3 || start < 1 || start >= end || end > p) {
        fprintf(
            stderr,
            "require prime>3 and 1<=start<end<=prime; got %d [%d,%d)\n",
            p, start, end
        );
        return 2;
    }

    Candidate *items = NULL;
    size_t count = 0, capacity = 0;
    uint64_t visited = 0;
    uint64_t categories[2][2] = {{0, 0}, {0, 0}};

    for (int p0 = start; p0 < end; ++p0) {
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

    printf("{\n");
    printf("  \"prime\": %d,\n", p);
    printf("  \"normalization\": \"e0=1 jet-sign representative\",\n");
    printf("  \"p0_start_inclusive\": %d,\n", start);
    printf("  \"p0_end_exclusive\": %d,\n", end);
    printf("  \"visited_parameter_tuples\": %llu,\n", (unsigned long long)visited);
    printf("  \"candidate_count\": %zu,\n", count);
    printf("  \"split_category_counts\": {\n");
    printf("    \"neither\": %llu,\n", (unsigned long long)categories[0][0]);
    printf("    \"iv_only\": %llu,\n", (unsigned long long)categories[0][1]);
    printf("    \"i4_only\": %llu,\n", (unsigned long long)categories[1][0]);
    printf("    \"both\": %llu\n", (unsigned long long)categories[1][1]);
    printf("  },\n");
    printf("  \"tuples\": [\n");
    for (size_t index = 0; index < count; ++index) {
        Candidate value = items[index];
        printf(
            "    [%d, %d, %d, %d, %d, %d, %d, %d]%s\n",
            value.p0, value.p1, value.p2, value.p3, value.r, value.s,
            value.split_i4, value.split_iv,
            index + 1 == count ? "" : ","
        );
    }
    printf("  ]\n");
    printf("}\n");

    free(items);
    return 0;
}
