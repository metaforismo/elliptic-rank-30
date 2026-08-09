#!/usr/bin/env python3
"""Derive a second rational split parameter by a tangent on the auxiliary quartic."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import sympy as sp

COEFFS=(172682352793305664,-537299076125816880,501757910513324641,-138642489541559748,-1600508800995452)
BASE_A,BASE_B,BASE_Y=34,19,79814494860
T_NUM=(-4386960832,7878955272,-3751492583,1932817518)
T_DEN=(658473248,6173506772,-12609842723,-1342798922)

def canonical_hash(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()

def build():
    s,u,X=sp.symbols('s u X'); A,B,C,D,E=map(sp.Integer,COEFFS)
    Q=A*s**4+B*s**3+C*s**2+D*s+E
    s0=sp.Rational(BASE_A,BASE_B); y0=sp.Rational(BASE_Y,BASE_B**2)
    assert sp.expand(Q.subs(s,s0)-y0**2)==0
    shifted=sp.Poly(sp.cancel((Q.subs(s,s0+u)-y0**2)/u),u,domain=sp.QQ)
    cubic=sp.expand(u*X**2+2*y0*X-shifted.as_expr())
    X0=sp.cancel(shifted.eval(0)/(2*y0))
    slope=sp.cancel(-sp.diff(cubic,u).subs({u:0,X:X0})/sp.diff(cubic,X).subs({u:0,X:X0}))
    tangent=sp.expand(X0+slope*u); intersection=sp.factor(cubic.subs(X,tangent))
    quotient=sp.Poly(sp.cancel(intersection/u**2),u,domain=sp.QQ); assert quotient.degree()==1
    u1=sp.cancel(-quotient.nth(0)/quotient.nth(1)); X1=sp.cancel(tangent.subs(u,u1)); s1=sp.cancel(s0+u1); y1=sp.cancel(y0+u1*X1)
    assert sp.expand(y1**2-Q.subs(s,s1))==0
    a1,b1=map(sp.Integer,sp.fraction(s1)); Y1=sp.cancel(y1*b1**2); assert Y1.is_Integer
    assert A*a1**4+B*a1**3*b1+C*a1**2*b1**2+D*a1*b1**3+E*b1**4==Y1**2
    def poly(c): return sum(sp.Integer(v)*s**(3-i) for i,v in enumerate(c))
    tn,td=poly(T_NUM),poly(T_DEN); t1=sp.cancel(tn.subs(s,s1)/td.subs(s,s1))
    cover=sp.Poly(sp.together(tn-t1*td),s,domain=sp.QQ); roots=[]
    for factor,multiplicity in sp.factor_list(cover.as_expr())[1]:
        p=sp.Poly(factor,s,domain=sp.QQ); assert p.degree()==1 and multiplicity==1; roots.append(sp.cancel(-p.nth(0)/p.nth(1)))
    roots=sorted(roots,key=lambda q:(int(sp.denom(q)),int(sp.numer(q))))
    assert len(roots)==3 and len(set(roots))==3 and s1 in roots and all(sp.cancel(tn.subs(s,r)/td.subs(s,r))==t1 for r in roots)
    payload={
      'schema_version':1,'certificate_id':'small-split-e8-three-trisection-tangent-split-parameter-v1','claim_status':'proved exact rational/integer arithmetic',
      'quartic_coefficients':list(COEFFS),'starting_point':{'a':str(BASE_A),'b':str(BASE_B),'Y':str(BASE_Y),'s':str(s0),'y':str(y0)},
      'shifted_cubic':str(cubic),'point_at_u_zero':{'u':'0','X':str(X0)},'tangent_slope':str(slope),'tangent_intersection_factorization':str(intersection),
      'new_affine_point':{'u':str(u1),'X':str(X1),'s':str(s1),'y':str(y1)},'new_homogeneous_point':{'a':str(a1),'b':str(b1),'Y':str(Y1)},
      'common_fibre_parameter':str(t1),'three_source_roots':[str(r) for r in roots],
      'exact_checks':{'starting_point_on_quartic':True,'tangent_has_double_intersection_at_start':True,'new_point_on_quartic':True,'homogeneous_quartic_identity':True,'three_distinct_rational_cover_roots':True,'all_roots_map_to_same_t':True},
      'truth_note':'This certifies a second completely split fibre of the common cubic cover. Mordell-Weil rank is a separate proof obligation.'}
    payload['record_sha256']=canonical_hash(payload); return payload

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--check',type=Path);args=ap.parse_args();p=build()
    if args.check:
        if json.loads(args.check.read_text())!=p:raise SystemExit('certificate mismatch')
        print('VERIFIED tangent split parameter sha256='+p['record_sha256'])
    else:
        text=json.dumps(p,indent=2,sort_keys=True)+'\n';args.output.write_text(text) if args.output else print(text,end='')
if __name__=='__main__':main()
