#!/usr/bin/env sage-python
"""Try exact Mordell-Weil rank computations after good finite-field reduction
for the invariant and anti-invariant rational elliptic surfaces in the marked
cyclic j=0 family.  A successful small fixed rank gives an arithmetic upper
bound in characteristic zero via Neron-Severi specialization.
"""
from sage.all import GF, FunctionField, EllipticCurve, QQ
import argparse, json
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('--prime',type=int,required=True)
ap.add_argument('--mu',default='3')
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args()
p=a.prime
F=GF(p)
K=FunctionField(F,'s'); s=K.gen()
mu=QQ(a.mu); q0=mu^2-mu+1; r0=mu*(mu-1)
aa=F(q0); bb=F(q0^3/r0)
B=bb^2*s^2-aa^3*(s+1)^3
curves={
  'invariant': EllipticCurve(K,[0,0,0,0,B]),
  'anti_invariant_standard_twist': EllipticCurve(K,[0,0,0,0,(4*s+1)^3*B]),
}
out={'prime':p,'mu':str(mu),'status':'pass','curves':{}}
for name,E in curves.items():
    row={'discriminant':str(E.discriminant())}
    try:
        row['rank']=int(E.rank())
    except Exception as exc:
        row['rank_error']=repr(exc)
    try:
        gs=E.gens()
        row['generator_count']=len(gs)
        row['generators']=[str(P) for P in gs]
    except Exception as exc:
        row['gens_error']=repr(exc)
    out['curves'][name]=row
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,sort_keys=True))
