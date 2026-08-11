/*
 * Exhaustive F_5 search for rational sections with P.O=2 on the five
 * translation-unique split semistable seed surfaces certified by the
 * rank17-split-seed-5 artifact.
 *
 * A section is represented as
 *
 *   x = X(t)/D(t)^2,  y = Y(t)/D(t)^3,
 *
 * with D monic quadratic, deg X<=8, deg Y<=12.  The search imposes the
 * nonidentity component conditions at I4 (t=0) and I3 (t=1), and rejects
 * common factors between D and X.  Every output is rechecked by the exact
 * polynomial identity over F_5.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define P 5
#define MAXD 24
#define MAX_RESULTS 4096

typedef struct {
    int c4[9];
    int c6[13];
    int source_parameters[8];
    int e0;
    int e1;
} Surface;

typedef struct {
    int surface_index;
    int D[3];
    int X[9];
    int Y[13];
    int component_tangent_sign_0;
    int component_tangent_sign_1;
} Section;

static const Surface surfaces[5] = {
    {{4,0,1,4,1,1,1,1,1},{2,0,2,3,0,3,3,0,0,0,0,4,1},{2,0,4,1,2,1,3,1},-1,1},
    {{4,4,1,0,2,4,0,3,1},{2,3,4,4,1,4,0,0,3,4,4,2,1},{2,1,0,0,2,2,0,3},-1,-1},
    {{4,3,3,1,0,4,3,0,1},{3,4,1,3,4,0,1,3,4,1,2,0,1},{2,2,1,3,2,2,0,0},1,-1},
    {{4,1,1,3,3,3,0,3,1},{2,2,4,2,0,4,1,0,1,0,4,2,1},{2,4,0,2,2,2,0,3},-1,1},
    {{4,1,0,2,0,3,3,0,1},{3,3,3,4,3,0,2,1,4,2,2,0,1},{2,4,1,1,2,2,0,0},1,1}
};

static int modp(int64_t x) { x%=P; if(x<0)x+=P; return (int)x; }
static int addp(int a,int b){return modp(a+b);} static int subp(int a,int b){return modp(a-b);} static int mulp(int a,int b){return modp((int64_t)a*b);}
static int powp(int a,int n){int r=1;while(n){if(n&1)r=mulp(r,a);a=mulp(a,a);n>>=1;}return r;}
static int invp(int a){if(!a){fprintf(stderr,"zero inverse\n");exit(4);}return powp(a,P-2);}

static void conv(const int *a,int da,const int *b,int db,int *out){
    memset(out,0,(da+db+1)*sizeof(int));
    for(int i=0;i<=da;++i)for(int j=0;j<=db;++j)out[i+j]=addp(out[i+j],mulp(a[i],b[j]));
}
static int degree(const int *a,int d){while(d>=0 && a[d]==0)--d;return d;}
static int eval1(const int *a,int d){int r=0;for(int i=0;i<=d;++i)r=addp(r,a[i]);return r;}
static int deriv1(const int *a,int d){int r=0;for(int i=1;i<=d;++i)r=addp(r,mulp(i%P,a[i]));return r;}

static int poly_gcd_degree(const int *aa,int da,const int *bb,int db){
    int a[MAXD+1]={0},b[MAXD+1]={0},r[MAXD+1]={0};
    memcpy(a,aa,(da+1)*sizeof(int));memcpy(b,bb,(db+1)*sizeof(int));
    da=degree(a,da);db=degree(b,db);
    while(db>=0){
        memcpy(r,a,sizeof(r));int dr=da;int il=invp(b[db]);
        while(dr>=db){int sh=dr-db;int q=mulp(r[dr],il);for(int i=0;i<=db;++i)r[i+sh]=subp(r[i+sh],mulp(q,b[i]));dr=degree(r,dr);}
        memcpy(a,b,sizeof(a));da=db;memcpy(b,r,sizeof(b));db=dr;
    }
    return da;
}

static int square_roots(int value,int roots[2]){
    int n=0;value=modp(value);for(int r=0;r<P;++r)if(mulp(r,r)==value){roots[n++]=r;if(n==2)break;}return n;
}

static int polynomial_square_roots(const int R[25], int roots[2][13]){
    int dr=degree(R,24);if(dr<0){return 0;}if(dr&1)return 0;int m=dr/2;if(m>12)return 0;
    int lead[2];int nlead=square_roots(R[dr],lead);int count=0;
    for(int lr=0;lr<nlead;++lr){
        int y[13]={0};y[m]=lead[lr];if(y[m]==0)continue;
        int inv2lead=invp(mulp(2,y[m]));
        int ok=1;
        for(int r=1;r<=m;++r){
            int k=2*m-r;int sum=0;
            for(int i=m-r+1;i<=m-1;++i){int j=k-i;if(j>=0 && j<=m)sum=addp(sum,mulp(y[i],y[j]));}
            y[m-r]=mulp(subp(R[k],sum),inv2lead);
        }
        int yy[25]={0};conv(y,m,y,m,yy);
        for(int i=0;i<=24;++i)if(yy[i]!=R[i]){ok=0;break;}
        if(ok){memcpy(roots[count++],y,sizeof(y));if(count==2)break;}
    }
    return count;
}

static void curve_rhs(const Surface *S,const int D[3],const int X[9],int R[25]){
    int D2[5]={0},D4[9]={0},D6[13]={0};int X2[17]={0},X3[25]={0};int temp1[17]={0},temp2[25]={0},temp3[25]={0};
    conv(D,2,D,2,D2);conv(D2,4,D2,4,D4);conv(D4,8,D2,4,D6);
    conv(X,8,X,8,X2);conv(X2,16,X,8,X3);
    conv(S->c4,8,X,8,temp1);conv(temp1,16,D4,8,temp2);conv(S->c6,12,D6,12,temp3);
    for(int i=0;i<=24;++i)R[i]=subp(subp(X3[i],mulp(3,temp2[i])),mulp(2,temp3[i]));
}

static int section_local_checks(const Surface *S,const int D[3],const int X[9],const int Y[13],int *sign0,int *sign1){
    int D2[5]={0};conv(D,2,D,2,D2);
    int xs0=modp(-S->e0*S->source_parameters[0]);
    int xs1=modp(-S->e1*S->source_parameters[4]);
    if(X[0]!=mulp(xs0,D2[0]))return 0;
    if(eval1(X,8)!=mulp(xs1,eval1(D2,4)))return 0;
    if(Y[0]!=0 || eval1(Y,12)!=0)return 0;
    int u0=subp(X[1],mulp(xs0,D2[1]));
    int u1=subp(deriv1(X,8),mulp(xs1,deriv1(D2,4)));
    if(!u0 || !u1)return 0;
    int d0=D[0],d1=eval1(D,2);if(!d0||!d1)return 0;
    int ratio0=mulp(Y[1],invp(mulp(d0,u0)));
    int ratio1=mulp(deriv1(Y,12),invp(mulp(d1,u1)));
    int tangent0=modp(-3*S->e0*S->source_parameters[0]);
    int tangent1=modp(-3*S->e1*S->source_parameters[4]);
    if(mulp(ratio0,ratio0)!=tangent0 || mulp(ratio1,ratio1)!=tangent1)return 0;
    *sign0=ratio0;*sign1=ratio1;return 1;
}

static void print_array(const int *a,int n){putchar('[');for(int i=0;i<n;++i){if(i)printf(", ");printf("%d",a[i]);}putchar(']');}

int main(void){
    Section results[MAX_RESULTS];int nresults=0;long long tested=0,square_rhs=0,local_pass=0;
    for(int si=0;si<5;++si){const Surface *S=&surfaces[si];
        for(int d0=1;d0<P;++d0)for(int d1=0;d1<P;++d1){int D[3]={d0,d1,1};if(eval1(D,2)==0)continue;int D2[5]={0};conv(D,2,D,2,D2);
        int xs0=modp(-S->e0*S->source_parameters[0]);int xs1=modp(-S->e1*S->source_parameters[4]);int x0=mulp(xs0,D2[0]);int target_sum=mulp(xs1,eval1(D2,4));
        for(int x2=0;x2<P;++x2)for(int x3=0;x3<P;++x3)for(int x4=0;x4<P;++x4)for(int x5=0;x5<P;++x5)for(int x6=0;x6<P;++x6)for(int x7=0;x7<P;++x7)for(int x8=0;x8<P;++x8){
            if(x8==P-1)continue;int X[9]={x0,0,x2,x3,x4,x5,x6,x7,x8};int sum_without_x1=0;for(int i=0;i<9;++i)sum_without_x1=addp(sum_without_x1,X[i]);X[1]=subp(target_sum,sum_without_x1);++tested;
            if(poly_gcd_degree(D,2,X,8)>0)continue;int R[25]={0};curve_rhs(S,D,X,R);int roots[2][13]={{0}};int nr=polynomial_square_roots(R,roots);if(!nr)continue;++square_rhs;
            for(int ri=0;ri<nr;++ri){int sg0=0,sg1=0;if(!section_local_checks(S,D,X,roots[ri],&sg0,&sg1))continue;++local_pass;if(nresults>=MAX_RESULTS){fprintf(stderr,"increase MAX_RESULTS\n");return 5;}Section *Q=&results[nresults++];Q->surface_index=si;memcpy(Q->D,D,sizeof(Q->D));memcpy(Q->X,X,sizeof(Q->X));memcpy(Q->Y,roots[ri],sizeof(Q->Y));Q->component_tangent_sign_0=sg0;Q->component_tangent_sign_1=sg1;}
        }}
    }
    printf("{\n  \"prime\": 5,\n  \"truth_status\": \"VERIFIED COMPUTATION\",\n  \"tested_reduced_ansatzes\": %lld,\n  \"square_rhs_ansatzes\": %lld,\n  \"section_count_up_to_y_sign\": %d,\n  \"sections\": [\n",tested,square_rhs,nresults/2);
    for(int i=0;i<nresults;++i){Section *Q=&results[i];printf("    {\"surface_index\": %d, \"D\": ",Q->surface_index);print_array(Q->D,3);printf(", \"X\": ");print_array(Q->X,9);printf(", \"Y\": ");print_array(Q->Y,13);printf(", \"tangent_ratio_0\": %d, \"tangent_ratio_1\": %d}%s\n",Q->component_tangent_sign_0,Q->component_tangent_sign_1,i+1==nresults?"":",");}
    printf("  ],\n  \"limitations\": [\"A finite-field section need not lift to characteristic zero.\", \"The search covers the stated degree bounds and local component conditions only.\"]\n}\n");
    return 0;
}
