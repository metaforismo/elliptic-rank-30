/*
 * Exhaustively filter the normalized rank-17 semistable fibre locus by the
 * arithmetic splitness conditions at the I4 and I3 fibres.
 *
 * This translation unit reuses the exact invariant construction and automatic
 * differentiation from rank17_semistable_fiber_search.c.  The original main
 * function is renamed before inclusion; all static helpers then remain visible
 * here without duplicating proof-critical formulas.
 */
#define main rank17_semistable_fiber_search_original_main
#include "rank17_semistable_fiber_search.c"
#undef main

#define MAX_SURFACES 512

typedef struct {
    int parameters[NVARS];
    int e0;
    int e1;
    int c4[9];
    int c6[13];
    int residual[6];
    int tangent_square_root_0;
    int tangent_square_root_1;
    int jacobian_rank;
} SplitSurface;

static int square_root_mod_or_minus_one(int value) {
    value = modp(value);
    for (int root = 0; root < P; ++root) {
        if (mulp(root, root) == value) return root;
    }
    return -1;
}

static int same_curve(
    const SplitSurface *surface,
    const int c4[9],
    const int c6[13]
) {
    return memcmp(surface->c4, c4, 9 * sizeof(int)) == 0
        && memcmp(surface->c6, c6, 13 * sizeof(int)) == 0;
}

static void print_surface(const SplitSurface *surface) {
    printf("      {\n");
    printf("        \"parameters\": ");
    print_int_array(surface->parameters, NVARS);
    printf(",\n");
    printf("        \"e0\": %d,\n", surface->e0 == P - 1 ? -1 : 1);
    printf("        \"e1\": %d,\n", surface->e1 == P - 1 ? -1 : 1);
    printf("        \"c4_coefficients_ascending\": ");
    print_int_array(surface->c4, 9);
    printf(",\n");
    printf("        \"c6_coefficients_ascending\": ");
    print_int_array(surface->c6, 13);
    printf(",\n");
    printf("        \"residual_quintic_coefficients_ascending\": ");
    print_int_array(surface->residual, 6);
    printf(",\n");
    printf("        \"split_tangent_square_root_at_0\": %d,\n",
           surface->tangent_square_root_0);
    printf("        \"split_tangent_square_root_at_1\": %d,\n",
           surface->tangent_square_root_1);
    printf("        \"formal_jacobian_rank\": %d\n",
           surface->jacobian_rank);
    printf("      }");
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

    long long raw_parameter_tuples = 0;
    long long exact_squarefree_tuples = 0;
    long long split_exact_squarefree_tuples = 0;
    long long split_by_sign[2][2] = {{0, 0}, {0, 0}};
    SplitSurface surfaces[MAX_SURFACES];
    int unique_surface_count = 0;

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
            ) continue;

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
                ) continue;

                int discriminant[25] = {0};
                if (!exact_fiber_orders(c4, c6, discriminant)) continue;
                int residual[6] = {0};
                int gcd_degree = residual_quintic_and_gcd_degree(
                    discriminant, residual
                );
                if (gcd_degree != 0) continue;
                ++exact_squarefree_tuples;

                Dual equations[NEQS];
                int derivative_c4[9], derivative_c6[13];
                build_equations(
                    values, e0, e1, equations,
                    derivative_c4, derivative_c6
                );
                int jacobian_rank = matrix_rank_6x8(equations);
                if (jacobian_rank != 6) {
                    fprintf(stderr, "unexpected nonsmooth exact point\n");
                    return 4;
                }

                /* At a signed multiplicative fibre the node is
                 * x=-e*l and the tangent cone is y^2=(-3*e*l)u^2.
                 * Rational nonidentity components therefore require
                 * -3*e*l to be a nonzero square in the residue field.
                 */
                int tangent0 = modp(-3LL * e0 * p0);
                int tangent1 = modp(-3LL * e1 * q0);
                int root0 = square_root_mod_or_minus_one(tangent0);
                int root1 = square_root_mod_or_minus_one(tangent1);
                if (root0 <= 0 || root1 <= 0) continue;
                ++split_exact_squarefree_tuples;
                ++split_by_sign[sign0][sign1];

                int duplicate = 0;
                for (int index = 0; index < unique_surface_count; ++index) {
                    if (same_curve(&surfaces[index], c4, c6)) {
                        duplicate = 1;
                        break;
                    }
                }
                if (duplicate) continue;
                if (unique_surface_count >= MAX_SURFACES) {
                    fprintf(stderr, "increase MAX_SURFACES\n");
                    return 5;
                }
                SplitSurface *surface = &surfaces[unique_surface_count++];
                memcpy(surface->parameters, values, sizeof(values));
                surface->e0 = e0;
                surface->e1 = e1;
                memcpy(surface->c4, c4, sizeof(c4));
                memcpy(surface->c6, c6, sizeof(c6));
                memcpy(surface->residual, residual, sizeof(residual));
                surface->tangent_square_root_0 = root0;
                surface->tangent_square_root_1 = root1;
                surface->jacobian_rank = jacobian_rank;
            }
        }
    }

    printf("{\n");
    printf("  \"schema_version\": 1,\n");
    printf("  \"truth_status\": \"VERIFIED COMPUTATION\",\n");
    printf("  \"prime\": %d,\n", P);
    printf("  \"raw_parameter_tuples\": %lld,\n", raw_parameter_tuples);
    printf("  \"smooth_squarefree_semistable_sign_parameter_tuples\": %lld,\n",
           exact_squarefree_tuples);
    printf("  \"split_smooth_squarefree_sign_parameter_tuples\": %lld,\n",
           split_exact_squarefree_tuples);
    printf("  \"split_counts_by_sign\": [\n");
    for (int sign0 = 0; sign0 < 2; ++sign0) {
        for (int sign1 = 0; sign1 < 2; ++sign1) {
            printf("    {\"e0\": %d, \"e1\": %d, \"count\": %lld}%s\n",
                   sign0 ? -1 : 1,
                   sign1 ? -1 : 1,
                   split_by_sign[sign0][sign1],
                   (sign0 == 1 && sign1 == 1) ? "" : ",");
        }
    }
    printf("  ],\n");
    printf("  \"unique_split_curves\": %d,\n", unique_surface_count);
    printf("  \"surfaces\": [\n");
    for (int index = 0; index < unique_surface_count; ++index) {
        print_surface(&surfaces[index]);
        printf("%s\n", index + 1 == unique_surface_count ? "" : ",");
    }
    printf("  ],\n");
    printf("  \"limitations\": [\n");
    printf("    \"Split fibres are necessary, not sufficient, for the required rational section.\",\n");
    printf("    \"Finite-field points are deformation seeds and do not imply characteristic-zero rational points.\"\n");
    printf("  ]\n");
    printf("}\n");
    return 0;
}
