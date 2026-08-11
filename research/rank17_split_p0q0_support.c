/* Count split smooth semistable seed tuples by (e0,e1,p0,q0). */
#define main rank17_semistable_fiber_search_original_main
#include "rank17_semistable_fiber_search.c"
#undef main

static int is_nonzero_square(int value) {
    value = modp(value);
    if (value == 0) return 0;
    for (int root = 1; root < P; ++root) {
        if (mulp(root, root) == value) return 1;
    }
    return 0;
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

    long long counts[2][2][P][P];
    memset(counts, 0, sizeof(counts));
    long long total = 0;

    for (int p0 = 1; p0 < P; ++p0)
    for (int p1 = 0; p1 < P; ++p1)
    for (int p2 = 0; p2 < P; ++p2)
    for (int p3 = 0; p3 < P; ++p3)
    for (int q0 = 1; q0 < P; ++q0)
    for (int q1 = 0; q1 < P; ++q1)
    for (int q2 = 0; q2 < P; ++q2)
    for (int s = 0; s < P; ++s) {
        int values[NVARS] = {p0,p1,p2,p3,q0,q1,q2,s};
        int c4[9], c6_base[13];
        build_c4_int(values, c4);
        build_c6_base_int(c4, c6_base);
        for (int sign0 = 0; sign0 < 2; ++sign0) {
            int e0 = sign0 ? P - 1 : 1;
            int target_low1 = mulp(e0, modp(3LL*p0*p0%P*p1));
            int target_low2 = mulp(e0, modp(3LL*p0*p0%P*p2 + 3LL*p0*p1%P*p1));
            int target_low3 = mulp(e0, modp(3LL*p0*p0%P*p3 + 6LL*p0*p1%P*p2 + (int64_t)p1*p1%P*p1));
            if (c6_base[1] != target_low1 || c6_base[2] != target_low2 || c6_base[3] != target_low3) continue;
            int c6[13];
            memcpy(c6, c6_base, sizeof(c6));
            c6[0] = mulp(e0, mulp(p0, mulp(p0,p0)));
            int value0=0,value1=0,value2=0;
            for (int i=0;i<=12;++i) {
                value0=addp(value0,c6[i]);
                value1=modp(value1+(int64_t)(i%P)*c6[i]);
                value2=modp(value2+(int64_t)binom_mod(i,2)*c6[i]);
            }
            for (int sign1=0;sign1<2;++sign1) {
                int e1=sign1?P-1:1;
                if (value0 != mulp(e1,mulp(q0,mulp(q0,q0)))) continue;
                if (value1 != mulp(e1,modp(3LL*q0*q0%P*q1))) continue;
                if (value2 != mulp(e1,modp(3LL*q0*q0%P*q2 + 3LL*q0*q1%P*q1))) continue;
                int discriminant[25]={0};
                if (!exact_fiber_orders(c4,c6,discriminant)) continue;
                int residual[6]={0};
                if (residual_quintic_and_gcd_degree(discriminant,residual) != 0) continue;
                Dual equations[NEQS]; int dc4[9],dc6[13];
                build_equations(values,e0,e1,equations,dc4,dc6);
                if (matrix_rank_6x8(equations) != 6) continue;
                if (!is_nonzero_square(modp(-3LL*e0*p0))) continue;
                if (!is_nonzero_square(modp(-3LL*e1*q0))) continue;
                ++counts[sign0][sign1][p0][q0];
                ++total;
            }
        }
    }

    printf("{\n  \"prime\": %d,\n  \"truth_status\": \"VERIFIED COMPUTATION\",\n",P);
    printf("  \"total_split_smooth_squarefree_tuples\": %lld,\n",total);
    printf("  \"support\": [\n");
    int first=1;
    for (int sign0=0;sign0<2;++sign0)
    for (int sign1=0;sign1<2;++sign1)
    for (int p0=1;p0<P;++p0)
    for (int q0=1;q0<P;++q0) {
        if (!counts[sign0][sign1][p0][q0]) continue;
        if (!first) printf(",\n");
        first=0;
        printf("    {\"e0\": %d, \"e1\": %d, \"p0\": %d, \"q0\": %d, \"count\": %lld}",
               sign0?-1:1, sign1?-1:1, p0,q0,counts[sign0][sign1][p0][q0]);
    }
    printf("\n  ]\n}\n");
    return 0;
}
