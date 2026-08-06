#!/usr/bin/env sage-python
"""Deterministically search rational cubic pencils with nine rational base
points and square-free degree-12 discriminant.

Such a pencil has twelve I1 fibres, no reducible fibres, and its blow-up is a
rational elliptic surface with split Mordell-Weil rank 8 over Q(t).  This is an
exact geometry certificate; no analytic rank heuristic is used.
"""
from sage.all import *
import argparse, json, random, traceback
from pathlib import Path


def cubic_basis(R):
    X,Y,Z=R.gens()
    return [X^3,X^2*Y,X^2*Z,X*Y^2,X*Y*Z,X*Z^2,Y^3,Y^2*Z,Y*Z^2,Z^3]


def point_sets(limit, seed):
    rng=random.Random(seed)
    anchor=[(0,0,1),(1,0,1),(0,1,1),(1,2,1)]
    pool=[(x,y,1) for x in range(-5,7) for y in range(-5,7)
          if (x,y,1) not in anchor]
    for _ in range(limit):
        yield anchor+rng.sample(pool,4)


def convert_cubic(cubic, O):
    errors=[]
    try:
        from sage.schemes.elliptic_curves.constructor import EllipticCurve_from_cubic
        return EllipticCurve_from_cubic(cubic,list(O)), 'EllipticCurve_from_cubic', errors
    except Exception as exc:
        errors.append('constructor:'+repr(exc))
    try:
        C=Curve(cubic)
        return C.elliptic_curve(list(O)), 'Curve.elliptic_curve', errors
    except Exception as exc:
        errors.append('curve_method:'+repr(exc))
    raise RuntimeError('; '.join(errors))


def normalise_poly(poly):
    poly=poly.univariate_polynomial()
    if not poly:
        return poly
    return poly/poly.lc()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--limit',type=int,default=60)
    ap.add_argument('--seed',type=int,default=20260807)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    R=PolynomialRing(QQ,names=('X','Y','Z')); X,Y,Z=R.gens()
    mons=cubic_basis(R)
    K=FunctionField(QQ,'t'); t=K.gen()
    RK=PolynomialRing(K,names=('X','Y','Z')); XX,YY,ZZ=RK.gens()
    attempts=[]
    winner=None
    for idx,pts in enumerate(point_sets(args.limit,args.seed),1):
        row={'attempt':idx,'points':[list(P) for P in pts]}
        try:
            M=Matrix(QQ,[[m(*P) for m in mons] for P in pts])
            row['evaluation_rank']=int(M.rank())
            if M.rank()!=8:
                row['status']='evaluation_rank_not_8'; attempts.append(row); continue
            ker=M.right_kernel()
            if ker.dimension()!=2:
                row['status']='not_a_pencil'; attempts.append(row); continue
            vv=ker.basis()
            F=sum(vv[0][i]*mons[i] for i in range(10))
            G=sum(vv[1][i]*mons[i] for i in range(10))

            # Find all rational affine intersections and the residual ninth point.
            A=PolynomialRing(QQ,names=('x','y'),order='lex'); x,y=A.gens()
            I=A.ideal([A(F(x,y,1)),A(G(x,y,1))])
            if I.dimension()!=0 or I.vector_space_dimension()!=9:
                row['status']='bad_intersection_scheme'
                row['intersection_dimension']=int(I.dimension())
                attempts.append(row); continue
            aff=[(QQ(v[x]),QQ(v[y]),QQ(1)) for v in I.variety(QQ)]
            residual=[P for P in aff if P not in set(pts)]
            if len(aff)!=9 or len(residual)!=1:
                row['status']='residual_not_unique_rational'
                row['rational_intersection_count']=len(aff)
                attempts.append(row); continue
            P9=residual[0]

            def lift(poly):
                return RK(poly(XX,YY,ZZ))
            cubic=lift(F)+t*lift(G)
            E,method,errors=convert_cubic(cubic,pts[0])
            delta=E.discriminant()
            num=delta.numerator(); den=delta.denominator()
            PR=num.parent()
            num=normalise_poly(PR(num))
            den=normalise_poly(PR(den))
            row.update({
              'F':str(F),'G':str(G),'ninth_point':[str(c) for c in P9],
              'conversion_method':method,'conversion_errors':errors,
              'a_invariants':[str(c) for c in E.a_invariants()],
              'discriminant_numerator':str(num),
              'discriminant_denominator':str(den),
              'discriminant_degree':int(num.degree()),
              'discriminant_squarefree':bool(num.gcd(num.derivative()).degree()==0),
              'discriminant_factorization':[(str(f),int(e)) for f,e in num.factor()],
            })
            if num.degree()==12 and row['discriminant_squarefree'] and den.degree()==0:
                row['status']='split_E8_certified'
                row['geometric_rank']=8
                row['arithmetic_rank_over_Qt']=8
                row['proof']='rational elliptic surface; 12 simple discriminant roots imply 12 I1 fibres; Shioda-Tate gives 10-2=8; rational base points split the section lattice'
                winner=row; attempts.append(row); break
            row['status']='not_twelve_I1'
        except Exception as exc:
            row['status']='error'; row['error']=repr(exc); row['traceback']=traceback.format_exc()
        attempts.append(row)
    out={'status':'pass' if winner else 'no_winner','seed':args.seed,'limit':args.limit,
         'attempts':attempts,'winner':winner}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':out['status'],'attempts':len(attempts),
                      'winner_attempt':winner and winner['attempt']},sort_keys=True))
    return 0 if winner else 1

if __name__=='__main__':
    raise SystemExit(main())
