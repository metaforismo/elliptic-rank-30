#!/usr/bin/env python3
"""Dependency-free exact certificate for the split E8 rational elliptic surface
in Kumar--Shioda, Example 18.

It verifies the eight polynomial sections, a square-free degree-11 finite
Discriminant (the twelfth I1 fibre is at infinity), and rank 8 of a specialized
subgroup at t=0 via exact finite-field mod-2 quotients.  Since specialization
cannot increase relations, the eight function-field sections are independent.
"""
from __future__ import annotations
from fractions import Fraction as Q
from functools import reduce
from math import gcd
import argparse, json
from pathlib import Path
from typing import Optional

mu=9699690
p2=Q(146156773903879871001810589,2**9*3*mu**2)
p1=-Q(24909805041567866985469379779685360019313,2**20*mu**3)
p0=Q(14921071761102637668643191215755039801471771138867387,2**23*3*mu**4)
q4=-Q(2243374456559366834339,2**5*mu**2)
q3=Q(430800343129403388346226518246078567,2**11*mu**3)
q2=Q(72555101947649011127391733034984158462573146409905769,2**22*3**2*mu**4)
q1=-Q(1288109930551729133820743237846836849158406377255698116491924530489,2**29*3*mu**5)
q0=Q(8827176793323619929427303381485459401911918837196838709750423283443360357992650203,2**42*3**3*mu**6)

# (g,a,b,v) where x=g*t^2+a*t+b and v is the multiplicative specialization.
XDATA=[
(Q(3),-Q(99950606190359,620780160),Q(4325327557647488120209649813,2642523476911718400),Q(3)),
(Q(5,4),-Q(153332163637781,1655413760),Q(5414114237697608646836821,5138596941004800),Q(5)),
(Q(7,9),-Q(203120672689603,2793510720),Q(6943164348569130636788638639,7927570430735155200),Q(7)),
(Q(11,25),-Q(8564057914757,147804800),Q(115126372233675800396600989,155442557465395200),Q(11)),
(Q(13,36),-Q(347479008951469,6385167360),Q(157133607680949617374030405417,221971972060584345600),Q(13)),
(Q(17,64),-Q(1327421017414859,26486620160),Q(5942419292933021418457517303,8901131711702630400),Q(17)),
(Q(19,81),-Q(489830985359431,10056638592),Q(46685137201743696441477454951,71348133876616396800),Q(19)),
(Q(120,169),-Q(30706596009257,440806080),Q(76164443074828743662165466409,55823308449760051200),Q(15,2)),
]

# Polynomials are coefficient lists in ascending order over Q.
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return a or [Q(0)]
def add(a,b):
    n=max(len(a),len(b)); return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)])
def neg(a): return [-x for x in a]
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)
def scale(a,c): return trim([c*x for x in a])
def powp(a,n):
    out=[Q(1)]; base=a
    while n:
        if n&1: out=mul(out,base)
        base=mul(base,base); n//=2
    return out
def deriv(a): return trim([Q(i)*a[i] for i in range(1,len(a))]) if len(a)>1 else [Q(0)]
def divmodp(a,b):
    a=trim(a); b=trim(b)
    if b==[0]: raise ZeroDivisionError
    q=[Q(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and a!=[0]:
        k=len(a)-len(b); c=a[-1]/b[-1]; q[k]=c
        a=sub(a,[Q(0)]*k+scale(b,c))
    return trim(q),trim(a)
def gcdp(a,b):
    a,b=trim(a),trim(b)
    while b!=[0]: _,r=divmodp(a,b); a,b=b,r
    return scale(a,1/a[-1]) if a!=[0] else [0]
def evalp(a,t):
    r=Q(0)
    for c in reversed(a): r=r*t+c
    return r

def qstr(x):
    return str(x.numerator) if x.denominator==1 else f"{x.numerator}/{x.denominator}"

# Curve: y^2=x^3+t^2*x^2+(p0+p1*t+p2*t^2)x+q0+...+q4*t^4+t^5.
A2=[0,0,1]
A4=[p0,p1,p2]
A6=[q0,q1,q2,q3,q4,1]


def reconstruct(g,a,b,v):
    u=(v+1)/(v-1); h=g*u
    assert h*h==g**3+g**2
    c=(3*a*g*g+2*a*g+1)/(2*h)
    d=(q4+g*p2+3*b*g*g+(2*b+3*a*a)*g+a*a-c*c)/(2*h)
    e=(q3+a*p2+g*p1+6*a*b*g+2*a*b+a**3-2*c*d)/(2*h)
    assert 2*c*e+d*d==q2+b*p2+a*p1+g*p0+3*b*b*g+b*b+3*a*a*b
    assert 2*d*e==q1+b*p1+a*p0+3*a*b*b
    assert e*e==q0+b*p0+b**3
    x=[b,a,g]; y=[e,d,c,h]
    rhs=add(add(powp(x,3),mul(A2,powp(x,2))),add(mul(A4,x),A6))
    assert sub(powp(y,2),rhs)==[0]
    return h,c,d,e

# Finite-field group law for y^2=x^3+a2*x^2+a4*x+a6.
Point=Optional[tuple[int,int]]
def fmod(x,p): return (x.numerator%p)*pow(x.denominator%p,-1,p)%p
def ecadd(P,R,p,a2,a4):
    if P is None:return R
    if R is None:return P
    x1,y1=P; x2,y2=R
    if x1==x2 and (y1+y2)%p==0:return None
    if P==R:
        if y1%p==0:return None
        lam=(3*x1*x1+2*a2*x1+a4)*pow(2*y1%p,-1,p)%p
    else: lam=(y2-y1)*pow((x2-x1)%p,-1,p)%p
    x3=(lam*lam-a2-x1-x2)%p
    return x3,(lam*(x1-x3)-y1)%p
def all_points(p,a2,a4,a6):
    pts=[None]
    for x in range(p):
        rhs=(x**3+a2*x*x+a4*x+a6)%p
        for y in range(p):
            if y*y%p==rhs: pts.append((x,y))
    return pts

def quotient_rows(p,a2,a4,a6,Ps):
    pts=all_points(p,a2,a4,a6); doubles={ecadd(P,P,p,a2,a4) for P in pts}
    rem=set(pts); cosets=[]
    while rem:
        rep=next(iter(rem)); cos={ecadd(rep,D,p,a2,a4) for D in doubles}
        cosets.append(cos); rem-=cos
    q=len(cosets)
    if q==1:return [],len(pts),q
    z=next(i for i,c in enumerate(cosets) if None in c)
    def ci(P):return next(i for i,c in enumerate(cosets) if P in c)
    def cadd(i,j):return ci(ecadd(next(iter(cosets[i])),next(iter(cosets[j])),p,a2,a4))
    basis=[]; span={z:0}
    for i in range(q):
        if i in span:continue
        k=len(basis); basis.append(i)
        for j,m in list(span.items()): span[cadd(j,i)]=m|(1<<k)
    rows=[[0]*len(Ps) for _ in basis]
    for j,P in enumerate(Ps):
        mask=span[ci(P)]
        for i in range(len(rows)):rows[i][j]=(mask>>i)&1
    return rows,len(pts),q

def rank2(rows):
    if not rows:return 0
    a=[sum((b&1)<<j for j,b in enumerate(r)) for r in rows]; rank=0
    for col in range(max(x.bit_length() for x in a)):
        piv=next((i for i in range(rank,len(a)) if (a[i]>>col)&1),None)
        if piv is None:continue
        a[rank],a[piv]=a[piv],a[rank]
        for i in range(len(a)):
            if i!=rank and (a[i]>>col)&1:a[i]^=a[rank]
        rank+=1
    return rank


def certificate():
    sections=[]
    for i,(g,a,b,v) in enumerate(XDATA,1):
        h,c,d,e=reconstruct(g,a,b,v)
        sections.append({'index':i,'v':qstr(v),'x':[qstr(g),qstr(a),qstr(b)],'y':[qstr(h),qstr(c),qstr(d),qstr(e)]})

    # Delta=-b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6, a1=a3=0.
    b2=scale(A2,4); b4=scale(A4,2); b6=scale(A6,4); b8=sub(scale(mul(A2,A6),4),powp(A4,2))
    delta=add(add(neg(mul(powp(b2,2),b8)),scale(powp(b4,3),-8)),add(scale(powp(b6,2),-27),scale(mul(mul(b2,b4),b6),9)))
    c4=sub(powp(b2,2),scale(b4,24))
    assert len(delta)-1==11
    assert len(gcdp(delta,deriv(delta)))-1==0
    assert len(gcdp(delta,c4))-1==0

    # Specialize at t=0 and obtain exact local certificate.
    t0=Q(0); sa2=evalp(A2,t0); sa4=evalp(A4,t0); sa6=evalp(A6,t0)
    spts=[]
    for sec in sections:
        g,a,b=map(Q,sec['x']); h,c,d,e=map(Q,sec['y'])
        x=b; y=e
        assert y*y==x**3+sa2*x*x+sa4*x+sa6
        spts.append((x,y))
    row_primes=[29,41,43,47,59,61,79,89]
    torsion_primes=[29,31,37]
    allrows=[]; local=[]; orders={}
    for p in sorted(set(row_primes+torsion_primes)):
        a2,a4,a6=map(lambda z:fmod(z,p),(sa2,sa4,sa6))
        disc=(a2*a2*a4*a4-4*a4**3-4*a2**3*a6-27*a6*a6+18*a2*a4*a6)%p
        assert disc!=0
        Ps=[(fmod(x,p),fmod(y,p)) for x,y in spts]
        rows,order,q=quotient_rows(p,a2,a4,a6,Ps); orders[p]=order
        before=rank2(allrows)
        if p in row_primes:allrows.extend(rows)
        local.append({'p':p,'order':order,'quotient_order':q,'rows':rows if p in row_primes else [],'rank_before':before,'rank_after':rank2(allrows)})
    assert reduce(gcd,[orders[p] for p in torsion_primes])==1
    assert rank2(allrows)==8
    return {
      'status':'pass','source':'Kumar--Shioda, Multiplicative excellent families, Example 18',
      'surface':{'equation':'y^2=x^3+t^2*x^2+(p0+p1*t+p2*t^2)*x+q0+q1*t+q2*t^2+q3*t^3+q4*t^4+t^5',
                 'mu':mu,'p':[qstr(p0),qstr(p1),qstr(p2)],'q':[qstr(q0),qstr(q1),qstr(q2),qstr(q3),qstr(q4)],
                 'finite_discriminant_degree':11,'finite_discriminant_squarefree':True,'gcd_delta_c4':1,
                 'fiber_configuration':'12 I1 (eleven finite roots and infinity)'},
      'sections':sections,
      'specialization_certificate':{'t':'0','torsion_primes':torsion_primes,'torsion_orders':[orders[p] for p in torsion_primes],
          'torsion_gcd':1,'row_primes':row_primes,'combined_mod2_rank':rank2(allrows),'local':local},
      'conclusion':'the eight displayed Q(t)-sections are independent; the rational elliptic surface has Mordell-Weil rank exactly 8 and lattice E8',
      'truth_note':'This is a rank-8 function-field base surface, not a rank-30 curve over Q.'}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); a=ap.parse_args()
    out=certificate(); text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text)
    print(json.dumps({'status':out['status'],'sections':len(out['sections']),'mod2_rank':out['specialization_certificate']['combined_mod2_rank']},sort_keys=True))
if __name__=='__main__':main()
