# An explicit split \(E_8\) rational elliptic surface over \(\mathbb Q(t)\)

## Purpose

The rank-30 packet construction requires a rational elliptic surface whose
Mordell--Weil group is already split over \(\mathbb Q(t)\) and has rank eight,
while spending no Néron--Severi rank on reducible fibres. Random cubic-pencil
conversion proved to be an unreliable computational representation. The
multiplicative excellent \(E_8\) family of Kumar--Shioda supplies the needed
object directly.

This note gives a dependency-free exact certificate for Example 18 of
Kumar--Shioda. It reconstructs the complete \(y\)-coordinates of the eight
published sections, verifies every section identity, proves that the finite
discriminant has eleven simple roots and no common root with \(c_4\), and proves
independence by specialization at \(t=0\) and exact finite-field reduction.

## The surface

Let \(\mu=9699690\). Put

\[
E:\quad y^2=x^3+t^2x^2+(p_0+p_1t+p_2t^2)x
       +q_0+q_1t+q_2t^2+q_3t^3+q_4t^4+t^5,
\]

where the exact coefficients are stored in
`certificates/kumar_shioda_split_e8_certificate.json` and constructed directly
by `research/kumar_shioda_split_e8_certificate.py`.

The eight sections have

\[
x_i(t)=g_i t^2+a_i t+b_i,
\qquad y_i(t)=h_i t^3+c_i t^2+d_i t+e_i.
\]

For each published quadruple \((g_i,a_i,b_i,v_i)\), the certificate reconstructs

\[
u_i=\frac{v_i+1}{v_i-1},\qquad h_i=g_i u_i,
\]

and then obtains \(c_i,d_i,e_i\) recursively from coefficient comparison. All
seven coefficient equations of \(y_i^2-x_i^3-t^2x_i^2-a_4x_i-a_6=0\) are
verified exactly over \(\mathbb Q\).

## Fibre configuration

For \(a_1=a_3=0\), the standard invariants give

\[
\Delta=-b_2^2b_8-8b_4^3-27b_6^2+9b_2b_4b_6.
\]

The exact certificate proves

\[
\deg_t\Delta=11,
\qquad \gcd(\Delta,\Delta')=1,
\qquad \gcd(\Delta,c_4)=1.
\]

Thus there are eleven finite nodal fibres. The missing degree of the minimal
discriminant is at infinity, so the rational elliptic surface has fibre
configuration

\[
12I_1.
\]

There is therefore no reducible-fibre root lattice. Shioda--Tate gives the
geometric Mordell--Weil rank

\[
10-2=8.
\]

## Exact arithmetic independence certificate

Specialize at \(t=0\). All eight sections specialize to rational points on a
nonsingular elliptic curve over \(\mathbb Q\).

Good reduction at \(29,31,37\) gives group orders

\[
36,\qquad39,\qquad37,
\]

whose gcd is one. Hence the specialized curve has trivial rational torsion.

The images of the eight points in the finite quotients

\[
E(\mathbb F_p)/2E(\mathbb F_p),
\quad p=29,41,43,47,59,61,79,89,
\]

produce an exact binary matrix of rank eight. Consequently the specialized
points are \(\mathbb Z\)-independent. Any relation among the generic sections
would specialize to a relation among these points, so the eight sections are
independent over \(\mathbb Q(t)\).

Combining this lower bound with Shioda--Tate proves

\[
\boxed{\operatorname{rank}E(\mathbb Q(t))=8.}
\]

The Mordell--Weil lattice is the split \(E_8\) lattice described by the source
construction.

## Consequence for the rank-30 search

This removes the first bottleneck in the quartic packet programme. We now have
an explicit, verified rank-eight base surface with twelve irreducible singular
fibres and all eight directions rational over \(\mathbb Q(t)\).

The next exact search is not another random cubic-pencil conversion. For
rational branch triples \(\{a,b,c\}\), compute the three twist ranks

\[
E^{(t-a)(t-b)},\qquad
E^{(t-a)(t-c)},\qquad
E^{(t-b)(t-c)}.
\]

A packet satisfying

\[
r_{ab}+r_{ac}+r_{bc}\ge22
\]

would force rank at least thirty after the corresponding rational biquadratic
base change. Finite-field scans are used only to reject weak branch triples or
to prioritize lifts; a characteristic-zero promotion requires explicit
\(\mathbb Q(t)\)-sections and a separate independence certificate.

## Truth status

This is an exact rank-eight function-field construction and search-enabling
intermediate theorem. It is not a rank-thirty elliptic curve over \(\mathbb Q\).
