#!/usr/bin/env sage-python
"""Construct a split rank-8 rational elliptic surface from eight rational
base points of a cubic pencil.  The ninth point is recovered exactly; the
generic cubic is converted to Weierstrass form and the base points are mapped
to sections.  This is a construction experiment, not a rank claim unless the
height matrix is certified non-singular.
"""
from sage.all import *
import json, traceback
from pathlib import Path

OUT=Path('search/results/cubic_pencil_e8.json')
OUT.parent.mkdir(parents=True,exist_ok=True)
result={'status':'started'}
try:
    Q=QQ
    R=PolynomialRing(Q,names=('X','Y','Z'))
    X,Y,Z=R.gens()
    mons=[X^3,X^2*Y,X^2*Z,X*Y^2,X*Y*Z,X*Z^2,Y^3,Y^2*Z,Y*Z^2,Z^3]
    pts=[
      (0,0,1),(1,0,1),(0,1,1),(1,1,1),
      (2,3,1),(3,2,1),(-1,2,1),(2,-1,1),
    ]
    M=Matrix(Q,[[m(*P) for m in mons] for P in pts])
    ker=M.right_kernel()
    result['evaluation_rank']=int(M.rank())
    result['pencil_dimension']=int(ker.dimension())
    if ker.dimension()!=2:
        raise RuntimeError('chosen points do not define a cubic pencil')
    vv=ker.basis()
    F=sum(vv[0][i]*mons[i] for i in range(10))
    G=sum(vv[1][i]*mons[i] for i in range(10))
    result['F']=str(F); result['G']=str(G)

    # Recover the residual ninth point in the affine chart Z=1.
    A=PolynomialRing(Q,names=('x','y'),order='lex'); x,y=A.gens()
    f=A(F(x,y,1)); g=A(G(x,y,1))
    I=A.ideal([f,g])
    result['intersection_dimension']=int(I.dimension())
    result['intersection_degree']=int(I.vector_space_dimension()) if I.dimension()==0 else None
    variety=I.variety(Q)
    aff=[(Q(v[x]),Q(v[y]),Q(1)) for v in variety]
    result['rational_affine_intersections']=[[str(c) for c in P] for P in aff]
    known=set(pts)
    residual=[P for P in aff if P not in known]
    if len(residual)!=1:
        raise RuntimeError(f'expected one rational residual point, got {residual}')
    P9=residual[0]
    allpts=pts+[P9]
    result['ninth_point']=[str(c) for c in P9]

    K=FunctionField(Q,'t'); t=K.gen()
    RK=PolynomialRing(K,names=('X','Y','Z')); XX,YY,ZZ=RK.gens()
    def lift(poly):
        return RK(str(poly))
    cubic=lift(F)+t*lift(G)

    from sage.schemes.elliptic_curves.constructor import EllipticCurve_from_cubic
    conversion=EllipticCurve_from_cubic(cubic,list(allpts[0]),morphism=True)
    result['conversion_type']=str(type(conversion))
    if isinstance(conversion,tuple):
        E=conversion[0]
        maps=conversion[1:]
    else:
        E=conversion; maps=[]
    result['weierstrass_a_invariants']=[str(a) for a in E.a_invariants()]
    result['discriminant']=str(E.discriminant())
    result['map_count']=len(maps)
    result['map_reprs']=[repr(m)[:1000] for m in maps]

    mapped=[]
    for P in allpts[1:]:
        image=None
        for m in maps:
            try:
                q=m(list(P))
                if q in E:
                    image=q; break
            except Exception:
                pass
        mapped.append(image)
    result['mapped_section_count']=sum(P is not None for P in mapped)
    result['mapped_sections']=[str(P) if P is not None else None for P in mapped]

    if all(P is not None for P in mapped):
        H=E.height_pairing_matrix(mapped)
        result['height_matrix']=[[str(x) for x in row] for row in H.rows()]
        result['height_determinant']=str(H.det())
        result['certified_rank']=int(H.rank())
        result['status']='pass' if H.rank()==8 else 'rank_below_8'
    else:
        result['status']='conversion_needs_map_api_fix'
except Exception as exc:
    result['status']='error'
    result['error']=repr(exc)
    result['traceback']=traceback.format_exc()
OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:result.get(k) for k in ('status','pencil_dimension','ninth_point','mapped_section_count','certified_rank','error')},sort_keys=True))
