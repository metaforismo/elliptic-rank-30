#!/usr/bin/env sage-python
"""Solve the complete minimal polynomial-section scheme for

    E_mu: y^2 = x^3 + b^2 r(t)^2 - a^3 q(t)^3,
    q=t^2-t+1, r=t(t-1),

where a=q(mu), b=q(mu)^3/r(mu).  The six zeros are the fibre of the
Legendre j-map above j(mu).  For generic mu this is a rational elliptic
surface with six type-II fibres and geometric Mordell-Weil lattice E8.

No coefficient box is used.  The seven coefficient equations for
x=A t^2+B t+C and y=D t^3+E t^2+F t+G are solved as a zero-dimensional
scheme over QQ.  Only exact rational solutions are promoted.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from sage.all import QQ, PolynomialRing


def as_fraction_string(x):
    x=QQ(x)
    return str(x.numerator()) if x.denominator()==1 else f"{x.numerator()}/{x.denominator()}"


def solve(mu):
    started=time.time()
    mu=QQ(mu)
    q0=mu^2-mu+1
    r0=mu*(mu-1)
    if r0==0:
        raise ValueError("mu must not be 0 or 1")
    a=q0
    b=q0^3/r0

    # Lex order with A last makes the terminal eliminant preferentially univariate in A.
    R=PolynomialRing(QQ,names=('G','F','E','D','C','B','A'),order='lex')
    G,F,E,D,C,B,A=R.gens()
    S=PolynomialRing(R,'t')
    t=S.gen()
    q=t^2-t+1
    r=t*(t-1)
    Bpoly=b^2*r^2-a^3*q^3
    x=A*t^2+B*t+C
    y=D*t^3+E*t^2+F*t+G
    diff=y^2-x^3-Bpoly
    eqs=[R(diff[i]) for i in range(7)]
    I=R.ideal(eqs)

    result={
        'status':'started',
        'mu':as_fraction_string(mu),
        'a':as_fraction_string(a),
        'b':as_fraction_string(b),
        'j_legendre':as_fraction_string(256*q0^3/r0^2),
        'ideal_dimension':int(I.dimension()),
        'equation_count':len(eqs),
    }
    if I.dimension()!=0:
        result['status']='non_zero_dimensional'
        result['elapsed_seconds']=time.time()-started
        return result

    gb=I.groebner_basis()
    result['groebner_basis_size']=len(gb)
    result['scheme_degree']=int(I.vector_space_dimension())

    univariate=[]
    for poly in gb:
        used=[v for v in R.gens() if poly.degree(v)>0]
        if len(used)==1:
            fac=[(str(f),int(e)) for f,e in poly.factor()]
            univariate.append({'variable':str(used[0]),'degree':int(poly.degree()),'factors':fac})
    result['univariate_groebner_polynomials']=univariate

    # Sage/Singular exact rational variety enumeration.  Every returned point is rechecked.
    solutions=[]
    try:
        variety=I.variety(QQ)
    except Exception as exc:
        result['variety_error']=repr(exc)
        variety=[]
    for sol in variety:
        vals={str(v):QQ(sol[v]) for v in R.gens()}
        if any(e.subs(sol)!=0 for e in eqs):
            raise AssertionError('variety returned a non-solution')
        solutions.append({k:as_fraction_string(v) for k,v in vals.items()})

    # Deduplicate signs: x is the section class, y and -y are one pair.
    xclasses={}
    for sol in solutions:
        key=(sol['A'],sol['B'],sol['C'])
        xclasses.setdefault(key,[]).append((sol['D'],sol['E'],sol['F'],sol['G']))
    result['rational_solution_count']=len(solutions)
    result['rational_section_pair_count']=len(xclasses)
    result['rational_sections']=[
        {'x_coefficients':list(key),'y_coefficients':list(vals[0]),'sign_count':len(vals)}
        for key,vals in sorted(xclasses.items())
    ]
    result['status']='pass'
    result['elapsed_seconds']=time.time()-started
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mu',default='3')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    out=solve(QQ(args.mu))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:out.get(k) for k in ('status','mu','scheme_degree','rational_section_pair_count','elapsed_seconds')},sort_keys=True))
    return 0 if out['status']=='pass' else 1

if __name__=='__main__':
    raise SystemExit(main())
