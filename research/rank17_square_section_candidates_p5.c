/* Extract every square RHS in all P.O=2 denominator charts over F_5. */
#define main rank17_split_section_search_finite_chart_main
#include "rank17_split_section_search_p5.c"
#undef main

#define MAX_SQUARE_CANDIDATES 128

typedef struct {
    int surface_index;
    int chart;
    int D[3];
    int X[9];
    int Y[13];
    int y0;
    int y1;
    int u0;
    int u1;
    int tangent_ratio_0;
    int tangent_ratio_1;
    int target_tangent_square_0;
    int target_tangent_square_1;
    int required_component_match;
} SquareCandidate;

static void scan_denominator(
    int surface_index,
    const int D[3],
    int chart,
    SquareCandidate candidates[MAX_SQUARE_CANDIDATES],
    int *candidate_count
) {
    const Surface *surface = &surfaces[surface_index];
    if (D[0] == 0 || eval1(D, 2) == 0) return;
    int D2[5] = {0};
    conv(D, 2, D, 2, D2);
    int singular_x0 = modp(-surface->e0 * surface->source_parameters[0]);
    int singular_x1 = modp(-surface->e1 * surface->source_parameters[4]);
    int x0 = mulp(singular_x0, D2[0]);
    int target_sum = mulp(singular_x1, eval1(D2, 4));

    for (int x2 = 0; x2 < P; ++x2)
    for (int x3 = 0; x3 < P; ++x3)
    for (int x4 = 0; x4 < P; ++x4)
    for (int x5 = 0; x5 < P; ++x5)
    for (int x6 = 0; x6 < P; ++x6)
    for (int x7 = 0; x7 < P; ++x7)
    for (int x8 = 0; x8 < P; ++x8) {
        if (D[2] != 0 && x8 == modp(-D[2] * D[2])) continue;
        int X[9] = {x0, 0, x2, x3, x4, x5, x6, x7, x8};
        int sum_without_x1 = 0;
        for (int index = 0; index < 9; ++index) {
            sum_without_x1 = addp(sum_without_x1, X[index]);
        }
        X[1] = subp(target_sum, sum_without_x1);
        if (poly_gcd_degree(D, 2, X, 8) > 0) continue;

        int right_hand_side[25] = {0};
        curve_rhs(surface, D, X, right_hand_side);
        int roots[2][13] = {{0}};
        int root_count = polynomial_square_roots(right_hand_side, roots);
        for (int root_index = 0; root_index < root_count; ++root_index) {
            if (*candidate_count >= MAX_SQUARE_CANDIDATES) {
                fprintf(stderr, "increase MAX_SQUARE_CANDIDATES\n");
                exit(5);
            }
            SquareCandidate *candidate = &candidates[(*candidate_count)++];
            memset(candidate, 0, sizeof(*candidate));
            candidate->surface_index = surface_index;
            candidate->chart = chart;
            memcpy(candidate->D, D, sizeof(candidate->D));
            memcpy(candidate->X, X, sizeof(candidate->X));
            memcpy(candidate->Y, roots[root_index], sizeof(candidate->Y));
            candidate->y0 = roots[root_index][0];
            candidate->y1 = eval1(roots[root_index], 12);
            candidate->u0 = subp(X[1], mulp(singular_x0, D2[1]));
            candidate->u1 = subp(
                deriv1(X, 8),
                mulp(singular_x1, deriv1(D2, 4))
            );
            candidate->target_tangent_square_0 = modp(
                -3 * surface->e0 * surface->source_parameters[0]
            );
            candidate->target_tangent_square_1 = modp(
                -3 * surface->e1 * surface->source_parameters[4]
            );
            int d0 = D[0];
            int d1 = eval1(D, 2);
            if (candidate->u0 && d0) {
                candidate->tangent_ratio_0 = mulp(
                    roots[root_index][1],
                    invp(mulp(d0, candidate->u0))
                );
            } else {
                candidate->tangent_ratio_0 = -1;
            }
            if (candidate->u1 && d1) {
                candidate->tangent_ratio_1 = mulp(
                    deriv1(roots[root_index], 12),
                    invp(mulp(d1, candidate->u1))
                );
            } else {
                candidate->tangent_ratio_1 = -1;
            }
            candidate->required_component_match =
                candidate->y0 == 0
                && candidate->y1 == 0
                && candidate->u0 != 0
                && candidate->u1 != 0
                && candidate->tangent_ratio_0 >= 0
                && candidate->tangent_ratio_1 >= 0
                && mulp(
                    candidate->tangent_ratio_0,
                    candidate->tangent_ratio_0
                ) == candidate->target_tangent_square_0
                && mulp(
                    candidate->tangent_ratio_1,
                    candidate->tangent_ratio_1
                ) == candidate->target_tangent_square_1;
        }
    }
}

int main(void) {
    SquareCandidate candidates[MAX_SQUARE_CANDIDATES];
    int candidate_count = 0;
    for (int surface_index = 0; surface_index < 5; ++surface_index) {
        int D0[3] = {1, 0, 0};
        scan_denominator(surface_index, D0, 0, candidates, &candidate_count);
        for (int d0 = 1; d0 < P; ++d0) {
            int D1[3] = {d0, 1, 0};
            scan_denominator(surface_index, D1, 1, candidates, &candidate_count);
        }
        for (int d0 = 1; d0 < P; ++d0)
        for (int d1 = 0; d1 < P; ++d1) {
            int D2[3] = {d0, d1, 1};
            scan_denominator(surface_index, D2, 2, candidates, &candidate_count);
        }
    }

    printf("{\n  \"prime\": 5,\n  \"truth_status\": \"VERIFIED COMPUTATION\",\n");
    printf("  \"square_candidate_count_including_y_sign\": %d,\n", candidate_count);
    printf("  \"candidates\": [\n");
    for (int index = 0; index < candidate_count; ++index) {
        SquareCandidate *candidate = &candidates[index];
        printf("    {\"surface_index\": %d, \"affine_degree_D\": %d, \"D\": ",
               candidate->surface_index, candidate->chart);
        print_array(candidate->D, 3);
        printf(", \"X\": ");
        print_array(candidate->X, 9);
        printf(", \"Y\": ");
        print_array(candidate->Y, 13);
        printf(", \"y_at_0\": %d, \"y_at_1\": %d, \"normal_displacement_0\": %d, \"normal_displacement_1\": %d, \"tangent_ratio_0\": %d, \"tangent_ratio_1\": %d, \"target_tangent_square_0\": %d, \"target_tangent_square_1\": %d, \"required_component_match\": %s}%s\n",
               candidate->y0,
               candidate->y1,
               candidate->u0,
               candidate->u1,
               candidate->tangent_ratio_0,
               candidate->tangent_ratio_1,
               candidate->target_tangent_square_0,
               candidate->target_tangent_square_1,
               candidate->required_component_match ? "true" : "false",
               index + 1 == candidate_count ? "" : ",");
    }
    printf("  ]\n}\n");
    return 0;
}
