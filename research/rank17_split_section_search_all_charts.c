/*
 * Complete the P.O=2 section search by including all affine charts of a
 * homogeneous quadratic denominator Z(T,U):
 *
 *   deg D=2: both intersections with O are finite;
 *   deg D=1: one intersection is at infinity;
 *   deg D=0: both intersections are at infinity.
 *
 * The original search covered only deg D=2.  This program reuses its exact
 * polynomial and local-component arithmetic and replaces only the search loop.
 */
#define main rank17_split_section_search_finite_chart_main
#include "rank17_split_section_search_p5.c"
#undef main

static void consider_denominator(
    int surface_index,
    const int D[3],
    int denominator_chart,
    Section results[MAX_RESULTS],
    int *result_count,
    long long tested_by_chart[3],
    long long square_rhs_by_chart[3],
    long long local_pass_by_chart[3]
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
        /* If D has full affine degree two, the section is finite at infinity.
         * Reject the singular point x=-1 of the I12 fibre there.  If deg D<2,
         * the section meets O at infinity and this affine-node test does not
         * apply. */
        if (D[2] != 0 && x8 == modp(-D[2] * D[2])) continue;

        int X[9] = {x0, 0, x2, x3, x4, x5, x6, x7, x8};
        int sum_without_x1 = 0;
        for (int index = 0; index < 9; ++index) {
            sum_without_x1 = addp(sum_without_x1, X[index]);
        }
        X[1] = subp(target_sum, sum_without_x1);
        ++tested_by_chart[denominator_chart];

        if (poly_gcd_degree(D, 2, X, 8) > 0) continue;
        int right_hand_side[25] = {0};
        curve_rhs(surface, D, X, right_hand_side);
        int roots[2][13] = {{0}};
        int root_count = polynomial_square_roots(right_hand_side, roots);
        if (!root_count) continue;
        ++square_rhs_by_chart[denominator_chart];

        for (int root_index = 0; root_index < root_count; ++root_index) {
            int tangent_sign_0 = 0;
            int tangent_sign_1 = 0;
            if (!section_local_checks(
                    surface,
                    D,
                    X,
                    roots[root_index],
                    &tangent_sign_0,
                    &tangent_sign_1
                )) {
                continue;
            }
            ++local_pass_by_chart[denominator_chart];
            if (*result_count >= MAX_RESULTS) {
                fprintf(stderr, "increase MAX_RESULTS\n");
                exit(5);
            }
            Section *section = &results[(*result_count)++];
            section->surface_index = surface_index;
            memcpy(section->D, D, sizeof(section->D));
            memcpy(section->X, X, sizeof(section->X));
            memcpy(section->Y, roots[root_index], sizeof(section->Y));
            section->component_tangent_sign_0 = tangent_sign_0;
            section->component_tangent_sign_1 = tangent_sign_1;
        }
    }
}

int main(int argc, char **argv) {
    int surface_count = (int)(sizeof(surfaces) / sizeof(surfaces[0]));
    int first_surface = 0;
    int last_surface = surface_count - 1;
    if (argc == 2) {
        first_surface = atoi(argv[1]);
        last_surface = first_surface;
        if (first_surface < 0 || first_surface >= surface_count) return 2;
    } else if (argc != 1) {
        return 2;
    }

    Section results[MAX_RESULTS];
    int result_count = 0;
    long long tested_by_chart[3] = {0, 0, 0};
    long long square_rhs_by_chart[3] = {0, 0, 0};
    long long local_pass_by_chart[3] = {0, 0, 0};

    for (int surface_index = first_surface;
         surface_index <= last_surface;
         ++surface_index) {
        /* deg D=0: Z(T,U)=U^2. */
        {
            int D[3] = {1, 0, 0};
            consider_denominator(
                surface_index,
                D,
                0,
                results,
                &result_count,
                tested_by_chart,
                square_rhs_by_chart,
                local_pass_by_chart
            );
        }

        /* deg D=1: Z(T,U)=T*U+d0*U^2. */
        for (int d0 = 1; d0 < P; ++d0) {
            int D[3] = {d0, 1, 0};
            consider_denominator(
                surface_index,
                D,
                1,
                results,
                &result_count,
                tested_by_chart,
                square_rhs_by_chart,
                local_pass_by_chart
            );
        }

        /* deg D=2: Z(T,U)=T^2+d1*T*U+d0*U^2. */
        for (int d0 = 1; d0 < P; ++d0)
        for (int d1 = 0; d1 < P; ++d1) {
            int D[3] = {d0, d1, 1};
            consider_denominator(
                surface_index,
                D,
                2,
                results,
                &result_count,
                tested_by_chart,
                square_rhs_by_chart,
                local_pass_by_chart
            );
        }
    }

    printf("{\n");
    printf("  \"prime\": %d,\n", P);
    printf("  \"truth_status\": \"VERIFIED COMPUTATION\",\n");
    printf("  \"surface_range\": [%d, %d],\n", first_surface, last_surface);
    printf("  \"denominator_charts\": [\n");
    for (int chart = 0; chart < 3; ++chart) {
        printf(
            "    {\"affine_degree_D\": %d, \"zero_intersection_multiplicity_at_infinity\": %d, \"tested_reduced_ansatzes\": %lld, \"square_rhs_ansatzes\": %lld, \"local_component_passes\": %lld}%s\n",
            chart,
            2 - chart,
            tested_by_chart[chart],
            square_rhs_by_chart[chart],
            local_pass_by_chart[chart],
            chart == 2 ? "" : ","
        );
    }
    printf("  ],\n");
    printf("  \"section_count_up_to_y_sign\": %d,\n", result_count / 2);
    printf("  \"sections\": [\n");
    for (int index = 0; index < result_count; ++index) {
        Section *section = &results[index];
        printf("    {\"surface_index\": %d, \"D\": ", section->surface_index);
        print_array(section->D, 3);
        printf(", \"X\": ");
        print_array(section->X, 9);
        printf(", \"Y\": ");
        print_array(section->Y, 13);
        printf(", \"tangent_ratio_0\": %d, \"tangent_ratio_1\": %d}%s\n",
               section->component_tangent_sign_0,
               section->component_tangent_sign_1,
               index + 1 == result_count ? "" : ",");
    }
    printf("  ],\n");
    printf("  \"limitations\": [\n");
    printf("    \"A finite-field section need not lift to characteristic zero.\",\n");
    printf("    \"The search covers every homogeneous quadratic denominator chart, the stated numerator degree bounds, and the required I4/I3 component conditions.\"\n");
    printf("  ]\n");
    printf("}\n");
    return 0;
}
