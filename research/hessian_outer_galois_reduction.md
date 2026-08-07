# Hessian reduction for cyclic cubic trisections

**Author:** Francesco Giannicola  
**Truth status:** new intermediate theorem; no rank-30 curve is claimed.

## Marked cubic-surface setting

Let

\[
X=\{F=0\}\subset\mathbf P(V)\cong\mathbf P^3
\]

be a smooth cubic surface over a field of characteristic different from two
and three. Let `ell` be a 3-secant line not contained in `X`. A hyperplane
`H` not containing `ell` meets it in one point

\[
q=H\cap\ell.
\]

The plane cubic

\[
C_H=X\cap H
\]

maps with degree three to the pencil of planes through `ell`; equivalently it
is projected from `q`.

The previous trisection theorem reduces every minimal norm-six trisection on
the split `E8` rational elliptic surface to precisely this geometry.

## The polar form

Write the cubic form as

\[
F(v)=T(v,v,v)
\]

for its symmetric trilinear polarization. For `q notin X`, put

\[
a=F(q)=T(q,q,q),
\qquad
\lambda_q(u)=T(q,q,u),
\]

and define the symmetric bilinear form

\[
\boxed{
B_q(u,v)
=aT(q,u,v)-\lambda_q(u)\lambda_q(v).
}
\]

It descends naturally to the three-dimensional quotient `V/<q>`.

## Hyperplane criterion

Let `H=<q,W>`, where `W` is a two-dimensional subspace of `V/<q>`. In
coordinates `(x,y,z)` on `H` with `q=(1:0:0)`, the restricted cubic has the
form

\[
F_H=a x^3+x^2L(y,z)+xQ(y,z)+C(y,z).
\]

The translation

\[
x=X-\frac{L}{3a}
\]

removes the `X^2` term. The remaining mixed quadratic is

\[
Q-\frac{L^2}{3a}.
\]

Therefore projection from `q` is geometrically cyclic exactly when

\[
\boxed{3aQ-L^2=0.}
\]

In invariant form this is

\[
\boxed{B_q|_W=0.}
\]

When it holds, the plane cubic becomes

\[
aX^3+G_3(y,z)=0,
\]

and the deck transformation is `X -> zeta_3 X`.

## Hessian theorem

Choose a complement to `q` and write the matrix of `T(q,-,-)` as

\[
M_q=
\begin{pmatrix}
 a&\lambda^t\\
 \lambda&Q_q
\end{pmatrix}.
\]

The matrix of `B_q` on `V/<q>` is

\[
aQ_q-\lambda\lambda^t.
\]

The block determinant identity gives

\[
\boxed{
\det B_q=a^2\det M_q.
}
\]

The Hessian matrix of `F` at `q` is `6M_q`. Hence

\[
\boxed{
\det B_q=0
\iff q\in\operatorname{Hess}(X).
}
\]

Over an algebraic closure, a symmetric form on a three-dimensional space has
a two-dimensional totally isotropic subspace exactly when it is singular.
Consequently:

\[
\boxed{
\begin{array}{c}
\text{A hyperplane }H\text{ gives a geometrically cyclic projection}\
\text{from }q=H\cap\ell
\end{array}
\iff
\begin{array}{c}
q\in\ell\cap\operatorname{Hess}(X),\ q\notin X,\\
H/\langle q\rangle\text{ is a maximal isotropic plane for }B_q.
\end{array}
}
\]

This is an equivalence over the algebraic closure. Rationality of `q`, the
isotropic plane, the hyperplane, and the cyclic cover is an additional descent
condition.

## Finiteness and exact cardinality bounds

The Hessian of a cubic surface is a quartic surface. If `ell` is not contained
in the Hessian, then

\[
\#(\ell\cap\operatorname{Hess}(X))\le4
\]

with multiplicity.

For a simple Hessian point:

* if `rank(B_q)=2`, there are exactly two maximal isotropic planes over the
  algebraic closure;
* if `rank(B_q)=1`, the radical is the unique maximal isotropic plane;
* if `rank(B_q)=0`, every plane is isotropic, a highly special positive-
  dimensional exception.

Thus a generic minimal trace class has at most

\[
\boxed{4\times2=8}
\]

geometrically cyclic hyperplane sections.

## Search consequence

The previous finite search still treated each trace class as a four-dimensional
linear system. The Hessian theorem collapses it to a quartic computation on a
line:

1. construct the marked cubic surface `(X,ell)` for the trace class;
2. restrict the Hessian quartic of `X` to `ell`;
3. factor one binary quartic;
4. discard roots lying on `X` or over bad elliptic fibres;
5. compute the one or two maximal isotropic planes of `B_q`;
6. recover the corresponding hyperplanes and plane cubics;
7. test smoothness, rational descent, and positive rank;
8. canonicalize the resulting cover across trace classes.

For one maximal ternary trace code, the worst-case generic candidate count is
reduced from a continuous `P3` search to at most

\[
240\times8=1920
\]

explicit hyperplane sections. The first decisive experiment is to construct
the marked cubic surface for one trace class and verify that the line is not
contained in the Hessian. The second is to find one cyclic cover shared by two
distinct trace classes. Eleven Hermitian-independent classes over one
positive-rank cover would force generic rank thirty.

The symbolic coefficient and block-determinant identities are checked by
`research/hessian_outer_galois_certificate.py`.
