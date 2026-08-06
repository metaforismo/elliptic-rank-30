#!/usr/bin/env sage-python
"""Search finite fields for exact 31-zero Riemann--Roch norm identities.

A seed consists of an elliptic curve y^2=x^3+a*x+b and 31 affine points with
distinct x-coordinates and sum zero.  Linear algebra in L(31O) then recovers

    f=A(x)+y*B(x),  deg A<=15, deg B=14,

whose norm A^2-(x^3+a*x+b)B^2 splits into the 31 corresponding linear
factors.  Centered integer lifting is scored only as a heuristic; no
characteristic-zero claim is made without exact rational factorization and
independence.
"""
from sage.all import *
import argparse, json, random, traceback
from pathlib import Path


def centered(c,p):
    n=int(c)
    return n if n<=p//2 else n-p


def poly_coeffs(poly,n):
    return [poly[i] if i<=poly.degree() else poly.base_ring()(0) for i in range(n)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--prime',type=int,required=True)
    ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--curves',type=int,default=120)
    ap.add_argument('--subset-attempts',type=int,default=300)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(); p=a.prime; F=GF(p); rng=random.Random(a.seed)
    Rx=PolynomialRing(F,'x'); x=Rx.gen()
    hits=[]; best_lift=None; errors=[]
    for curve_index in range(a.curves):
        aa=F(rng.randrange(p)); bb=F(rng.randrange(p))
        E=EllipticCurve(F,[0,0,0,aa,bb])
        if E.discriminant()==0: continue
        affine=[P for P in E if not P.is_zero()]
        if len({int(P[0]) for P in affine})<31: continue
        for attempt in range(a.subset_attempts):
            rng.shuffle(affine)
            chosen=[]; used=set()
            for P in affine:
                xx=int(P[0])
                if xx in used: continue
                chosen.append(P); used.add(xx)
                if len(chosen)==30: break
            if len(chosen)<30: continue
            S=sum(chosen,E(0)); target=-S
            if target.is_zero() or int(target[0]) in used: continue
            points=chosen+[target]
            if len({int(P[0]) for P in points})!=31: continue
            assert sum(points,E(0)).is_zero()
            M=Matrix(F,[[P[0]^i for i in range(16)]+[P[1]*P[0]^i for i in range(15)] for P in points])
            ker=M.right_kernel()
            if ker.dimension()!=1: continue
            v=ker.basis()[0]
            if v[30]==0: continue
            v=v/v[30]
            A=sum(v[i]*x^i for i in range(16))
            B=sum(v[16+i]*x^i for i in range(15))
            R=x^3+aa*x+bb
            N=A^2-R*B^2
            if N.degree()!=31 or N.gcd(N.derivative()).degree()!=0: continue
            fac=N.factor()
            if any(f.degree()!=1 or e!=1 for f,e in fac): continue
            roots={int(-f[0]/f[1]) for f,e in fac}
            if roots!={int(P[0]) for P in points}: continue

            # Centered integer lift: exact over Q, scored by actual Q-linear factors.
            RQ=PolynomialRing(QQ,'x'); X=RQ.gen()
            ai=centered(aa,p); bi=centered(bb,p)
            AQ=sum(centered(v[i],p)*X^i for i in range(16))
            BQ=sum(centered(v[16+i],p)*X^i for i in range(15))
            NQ=AQ^2-(X^3+ai*X+bi)*BQ^2
            fQ=NQ.factor()
            linear=[]
            for ff,ee in fQ:
                if ff.degree()==1:
                    rr=-ff[0]/ff[1]
                    linear.extend([str(rr)]*ee)
            rec={
              'curve_index':curve_index,'subset_attempt':attempt,
              'a':int(aa),'b':int(bb),
              'A_coefficients':[int(c) for c in poly_coeffs(A,16)],
              'B_coefficients':[int(c) for c in poly_coeffs(B,15)],
              'x_roots':sorted(roots),
              'norm_leading_coefficient':int(N.leading_coefficient()),
              'finite_field_identity_verified':True,
              'centered_lift':{
                'a':ai,'b':bi,
                'A_coefficients':[centered(v[i],p) for i in range(16)],
                'B_coefficients':[centered(v[16+i],p) for i in range(15)],
                'rational_linear_factor_count':len(linear),
                'rational_linear_roots':linear,
                'factor_degrees':[(int(ff.degree()),int(ee)) for ff,ee in fQ],
              },
            }
            hits.append(rec)
            if best_lift is None or len(linear)>best_lift['centered_lift']['rational_linear_factor_count']:
                best_lift=rec
            break
    out={
      'status':'pass','prime':p,'seed':a.seed,'curves_tested':a.curves,
      'seed_count':len(hits),'best_centered_lift':best_lift,'seeds':hits,
      'truth_note':'finite-field norm identities are exact; centered Q lifts are heuristic until their norms split and points are certified',
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'prime':p,'seed_count':len(hits),
                      'best_Q_linear_factors':best_lift and best_lift['centered_lift']['rational_linear_factor_count']},sort_keys=True))

if __name__=='__main__': main()
