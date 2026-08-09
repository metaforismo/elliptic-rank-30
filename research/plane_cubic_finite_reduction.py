#!/usr/bin/env python3
"""Exact finite-reduction certificates directly on a smooth plane cubic.

The identity may be any rational point O, not necessarily a flex.  If
``third(P,Q)`` is the third intersection of the line PQ with the cubic
(tangent when P=Q), then

    P + Q = third(third(P,Q), O)

is the elliptic-curve group law with identity O.  This avoids a potentially
huge Weierstrass transformation when certifying specialized plane cubics.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence

Point=tuple[int,int,int]
Line=tuple[int,int,int]

F_COEFFS={
 (3,0,0):3972,(2,1,0):8080,(2,0,1):-65622,(1,2,0):31679,
 (1,1,1):-104467,(1,0,2):-232614,(0,3,0):24484,
 (0,2,1):-15556,(0,1,2):-173688,
}
G_COEFFS={
 (3,0,0):33084,(2,1,0):44912,(2,0,1):-62778,(1,2,0):24409,
 (1,1,1):-70613,(1,0,2):-138714,(0,3,0):-36220,
 (0,2,1):-122924,(0,1,2):347376,(0,0,3):1042128,
}

def normalize(values:Sequence[int],p:int)->tuple[int,...]:
    """Canonical projective normalization compatible with enumeration.

    Use the last nonzero coordinate as pivot: z=1 in the affine chart,
    otherwise y=1 at infinity, otherwise [1:0:0].
    """
    v=[int(x)%p for x in values]
    pivot=next((v[i] for i in range(len(v)-1,-1,-1) if v[i]),None)
    if pivot is None:raise ValueError('zero projective vector')
    z=pow(pivot,-1,p)
    return tuple((x*z)%p for x in v)

def line_through(P:Point,Q:Point,p:int)->Line:
    a=(P[1]*Q[2]-P[2]*Q[1])%p
    b=(P[2]*Q[0]-P[0]*Q[2])%p
    c=(P[0]*Q[1]-P[1]*Q[0])%p
    return normalize((a,b,c),p)  # type: ignore[return-value]

def line_value(L:Line,P:Point,p:int)->int:
    return sum(a*b for a,b in zip(L,P))%p

def combine_coeffs(t_num:int,t_den:int,p:int)->dict[tuple[int,int,int],int]:
    keys=set(F_COEFFS)|set(G_COEFFS)
    return {m:(t_den*F_COEFFS.get(m,0)+t_num*G_COEFFS.get(m,0))%p for m in keys}

def poly_value(coeffs:dict[tuple[int,int,int],int],P:Point,p:int)->int:
    x,y,z=P
    return sum(c*pow(x,i,p)*pow(y,j,p)*pow(z,k,p) for (i,j,k),c in coeffs.items())%p

def gradient(coeffs:dict[tuple[int,int,int],int],P:Point,p:int)->Line:
    x,y,z=P; out=[]
    for axis in range(3):
        value=0
        for exps,c in coeffs.items():
            e=exps[axis]
            if not e:continue
            powers=[x,y,z];term=c*e
            for q,a in enumerate(exps):
                term*=pow(powers[q],a-(1 if q==axis else 0),p)
            value+=term
        out.append(value%p)
    return normalize(out,p)  # type: ignore[return-value]

def enumerate_projective_points(p:int)->list[Point]:
    points=[(x,y,1) for x in range(p) for y in range(p)]
    points.extend((x,1,0) for x in range(p))
    points.append((1,0,0))
    return points

def rational_projective_mod(P:Sequence[Fraction|int],p:int)->Point:
    values=[]
    for q0 in P:
        q=Fraction(q0)
        if q.denominator%p==0:raise ZeroDivisionError('point denominator')
        values.append((q.numerator%p)*pow(q.denominator%p,-1,p)%p)
    return normalize(values,p)  # type: ignore[return-value]

@dataclass
class PlaneCubicGroup:
    p:int
    coeffs:dict[tuple[int,int,int],int]
    origin:Point
    points:list[Point]
    point_set:set[Point]
    tangents:dict[Point,Line]

    @classmethod
    def create(cls,p:int,t_num:int,t_den:int,origin_rational:Sequence[Fraction|int]):
        coeffs=combine_coeffs(t_num,t_den,p)
        points=[P for P in enumerate_projective_points(p) if poly_value(coeffs,P,p)==0]
        point_set=set(points)
        tangents={P:gradient(coeffs,P,p) for P in points}
        if any(all(v%p==0 for v in cls.raw_gradient(coeffs,P,p)) for P in points):
            raise ValueError('singular reduction')
        origin=rational_projective_mod(origin_rational,p)
        if origin not in point_set:raise ValueError('origin not on reduced cubic')
        group=cls(p,coeffs,origin,points,point_set,tangents)
        group.sanity_check()
        return group

    @staticmethod
    def raw_gradient(coeffs,P,p):
        x,y,z=P;out=[]
        for axis in range(3):
            value=0
            for exps,c in coeffs.items():
                e=exps[axis]
                if not e:continue
                powers=[x,y,z];term=c*e
                for q,a in enumerate(exps):term*=pow(powers[q],a-(1 if q==axis else 0),p)
                value+=term
            out.append(value%p)
        return tuple(out)

    def _add_vectors(self,P:Point,Q:Point,sign:int=1)->Point:
        return tuple((a+sign*b)%self.p for a,b in zip(P,Q))  # type: ignore[return-value]

    def _tangent_companion(self,P:Point)->Point:
        """Choose a second vector on the tangent line, independent of P."""
        a,b,c=self.tangents[P];p=self.p
        candidates=[]
        if a or b:candidates.append((b%p,(-a)%p,0))
        if a or c:candidates.append((c%p,0,(-a)%p))
        if b or c:candidates.append((0,c%p,(-b)%p))
        for V in candidates:
            if V!=(0,0,0) and normalize(V,p)!=P:
                return V
        raise AssertionError(('failed to choose tangent companion',P,self.tangents[P]))

    def third(self,P:Point,Q:Point)->Point:
        """Return the third intersection, with multiplicity, in O(1) field work.

        For P != Q, write the restricted binary cubic as
        A*s^2*t + B*s*t^2.  Its third root is [s:t]=[B:-A].
        For a tangent at P, choose V on the tangent and write it as
        A*s*t^2 + B*t^3; the third root is [-B:A].
        """
        p=self.p; inv2=pow(2,-1,p)
        if P!=Q:
            plus=poly_value(self.coeffs,self._add_vectors(P,Q,1),p)
            minus=poly_value(self.coeffs,self._add_vectors(P,Q,-1),p)
            A=((plus-minus)*inv2)%p
            B=((plus+minus)*inv2)%p
            if A==0 and B==0:
                raise AssertionError(('line contained in cubic',P,Q))
            return normalize(tuple((B*x-A*y)%p for x,y in zip(P,Q)),p)  # type: ignore[return-value]

        V=self._tangent_companion(P)
        plus=poly_value(self.coeffs,self._add_vectors(P,V,1),p)
        minus=poly_value(self.coeffs,self._add_vectors(P,V,-1),p)
        A=((plus+minus)*inv2)%p
        B=((plus-minus)*inv2)%p
        if B!=poly_value(self.coeffs,V,p):
            raise AssertionError(('tangent restriction mismatch',P,V,A,B))
        if A==0:
            return P
        return normalize(tuple((-B*x+A*y)%p for x,y in zip(P,V)),p)  # type: ignore[return-value]

    def add(self,P:Point,Q:Point)->Point:
        return self.third(self.third(P,Q),self.origin)

    def neg(self,P:Point)->Point:
        return self.third(self.third(self.origin,self.origin),P)

    def mul(self,n:int,P:Point)->Point:
        if n<0:return self.mul(-n,self.neg(P))
        R=self.origin;Q=P
        while n:
            if n&1:R=self.add(R,Q)
            n//=2
            if n:Q=self.add(Q,Q)
        return R

    def sanity_check(self):
        O=self.origin
        for P in self.points:
            if self.add(P,O)!=P or self.add(O,P)!=P:raise AssertionError('identity failure')
            if self.add(P,self.neg(P))!=O:raise AssertionError('inverse failure')
        n=len(self.points)
        for k in range(min(300,n*n)):
            P=self.points[(17*k+1)%n];Q=self.points[(31*k+3)%n];R=self.points[(47*k+5)%n]
            if self.add(self.add(P,Q),R)!=self.add(P,self.add(Q,R)):
                raise AssertionError('associativity sample failed')

    def quotient_vectors_points(self,reduced_points:Sequence[Point],ell:int):
        elements=self.points; H={self.mul(ell,P) for P in elements}; un=set(elements); cosets=[]; which={}
        while un:
            a=next(iter(un)); C={self.add(a,h) for h in H}; idx=len(cosets)
            for P in C:which[P]=idx
            un-=C;cosets.append(C)
        zero=which[self.origin]; reps=[next(iter(C)) for C in cosets]; span={zero:()}
        while len(span)<len(cosets):
            b=next(i for i in range(len(cosets)) if i not in span);new={}
            for i,v in span.items():
                for c in range(ell):new[which[self.add(reps[i],self.mul(c,reps[b]))]]=v+(c,)
            span=new
        for P in reduced_points:
            if P not in self.point_set:raise ValueError('reduced point off cubic')
        return [span[which[P]] for P in reduced_points],len(next(iter(span.values())))

    def quotient_vectors(self,rational_points:Sequence[Sequence[Fraction|int]],ell:int):
        reduced=[]
        for P0 in rational_points:
            P=rational_projective_mod(P0,self.p)
            if P not in self.point_set:raise ValueError('reduced point off cubic')
            reduced.append(P)
        return self.quotient_vectors_points(reduced,ell)

def row_rank(rows:Sequence[Sequence[int]],p:int)->int:
    if not rows:return 0
    A=[[x%p for x in row] for row in rows];m=len(A);n=len(A[0]);r=0
    for c in range(n):
        q=next((i for i in range(r,m) if A[i][c]),None)
        if q is None:continue
        A[r],A[q]=A[q],A[r];z=pow(A[r][c],-1,p);A[r]=[(z*x)%p for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                z=A[i][c];A[i]=[(x-z*y)%p for x,y in zip(A[i],A[r])]
        r+=1
        if r==m:break
    return r

def finite_reduction_certificate(*,t_num:int,t_den:int,origin:Sequence[Fraction|int],points:Sequence[Sequence[Fraction|int]],ell:int,primes:Iterable[int]):
    rows=[[] for _ in points];records=[];torsion=None
    for p in primes:
        if p<5:continue
        try:
            group=PlaneCubicGroup.create(p,t_num,t_den,origin)
            vectors,dim=group.quotient_vectors(points,ell)
        except (ValueError,ZeroDivisionError,AssertionError):
            continue
        order=len(group.points)
        if torsion is None and order%ell:torsion={'prime':p,'group_order':order}
        if dim:
            records.append({'prime':p,'group_order':order,'quotient_dimension':dim})
            for i,v in enumerate(vectors):rows[i].extend(v)
        rank=row_rank(rows,ell)
        if rank==len(points) and torsion is not None:break
    return {'ell':ell,'rank':row_rank(rows,ell),'rows':rows,'primes':records,'torsion_witness':torsion}
