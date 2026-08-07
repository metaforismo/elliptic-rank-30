#!/usr/bin/env sage-python
"""Finite-field laboratory for the rank-30 packet mechanism.

For deterministic cubic pencils through eight F_p-rational points:
  * convert the generic cubic to an elliptic curve over F_p(t),
  * require exact function-field rank 8 for the base surface,
  * compute exact ranks of every quadratic twist ramified at a pair from a
    three-point branch support,
  * retain the best V4 packet by 8+r_ab+r_ac+r_bc.

Finite-field rank is discovery evidence only; lifting and explicit sections are
required before any characteristic-zero claim.
"""
from sage.all import *
import argparse, itertools, json, random, traceback
from pathlib import Path


def cubic_basis(R):
    X,Y,Z=R.gens()
    return [X**3,X**2*Y,X**2*Z,X*Y**2,X*Y*Z,X*Z**2,Y**3,Y**2*Z,Y*Z**2,Z**3]


def convert(cubic,O):
    from sage.schemes.elliptic_curves.constructor import EllipticCurve_from_cubic
    return EllipticCurve_from_cubic(cubic,list(O))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--prime',type=int,required=True)
    ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--surfaces',type=int,default=12)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    F=GF(a.prime); rng=random.Random(a.seed)
    R=PolynomialRing(F,names=('X','Y','Z')); X,Y,Z=R.gens(); mons=cubic_basis(R)
    K=FunctionField(F,'t'); t=K.gen()
    RK=PolynomialRing(K,names=('X','Y','Z')); XX,YY,ZZ=RK.gens()
    affine=[(F(x),F(y),F(1)) for x in range(a.prime) for y in range(a.prime)]
    results=[]; best=None
    for sid in range(a.surfaces):
        row={'surface_index':sid}
        try:
            pts=rng.sample(affine,8)
            M=Matrix(F,[[m(*P) for m in mons] for P in pts])
            if M.rank()!=8:
                row['status']='bad_evaluation_rank'; results.append(row); continue
            ker=M.right_kernel()
            if ker.dimension()!=2:
                row['status']='not_pencil'; results.append(row); continue
            vv=ker.basis()
            F0=sum(vv[0][i]*mons[i] for i in range(10))
            G0=sum(vv[1][i]*mons[i] for i in range(10))
            cubic=RK(F0(XX,YY,ZZ))+t*RK(G0(XX,YY,ZZ))
            E=convert(cubic,pts[0])
            base_rank=int(E.rank())
            row['base_rank']=base_rank
            row['points']=[[int(c) for c in P] for P in pts]
            row['F']=str(F0); row['G']=str(G0)
            row['a_invariants']=[str(c) for c in E.a_invariants()]
            if base_rank!=8:
                row['status']='base_rank_not_8'; results.append(row); continue

            supports=list(F)+['infinity']
            pair_rank={}
            for aa,bb in itertools.combinations(supports,2):
                if aa=='infinity': aa,bb=bb,aa
                if bb=='infinity':
                    d=t-F(aa)
                else:
                    d=(t-F(aa))*(t-F(bb))
                try:
                    Et=E.quadratic_twist(d)
                    rr=int(Et.rank())
                except Exception as exc:
                    pair_rank[str((aa,bb))]={'error':repr(exc)}
                    continue
                pair_rank[str((aa,bb))]={'rank':rr,'d':str(d)}
            row['pair_ranks']=pair_rank

            packets=[]
            for triple in itertools.combinations(supports,3):
                pairs=list(itertools.combinations(triple,2))
                vals=[]; ok=True
                for pp in pairs:
                    x,y=pp
                    if x=='infinity': x,y=y,x
                    key=str((x,y))
                    rec=pair_rank.get(key,{})
                    if 'rank' not in rec:
                        ok=False; break
                    vals.append(rec['rank'])
                if ok:
                    total=8+sum(vals)
                    packets.append({'support':[str(x) for x in triple],
                                    'twist_ranks':vals,'total_rank':total})
            packets.sort(key=lambda z:(z['total_rank'],z['twist_ranks']),reverse=True)
            row['best_packets']=packets[:10]
            row['status']='pass'
            if packets and (best is None or packets[0]['total_rank']>best['packet']['total_rank']):
                best={'surface':row,'packet':packets[0]}
        except Exception as exc:
            row['status']='error'; row['error']=repr(exc); row['traceback']=traceback.format_exc()
        results.append(row)
    out={'status':'pass','prime':a.prime,'seed':a.seed,'surface_count':a.surfaces,
         'best':best,'results':results,
         'truth_note':'finite-field discovery evidence only; no Q(t) or rank-30 claim'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'prime':a.prime,'seed':a.seed,
                      'best_total':best and best['packet']['total_rank']},sort_keys=True))

if __name__=='__main__': main()
