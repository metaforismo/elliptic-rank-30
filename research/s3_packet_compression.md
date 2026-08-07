# S3 compression of a genus-zero rank packet

## The invariant family

Put

\[
q(t)=t^2-t+1,
\qquad
r(t)=t(t-1),
\]

and consider

\[
E:\quad y^2=x^3+B(t),
\qquad
B(t)=b^2r(t)^2-a^3q(t)^3.
\]

The transformations

\[
s(t)=1-t,
\qquad i(t)=1/t
\]

generate the anharmonic group \(S_3\) permuting \(0,1,\infty\). Exact identities give

\[
q(1-t)=q(t),\quad r(1-t)=r(t),\quad
t^2q(1/t)=q(t),\quad t^3r(1/t)=-r(t).
\]

Consequently

\[
B(1-t)=B(t),
\qquad t^6B(1/t)=B(t).
\]

The first identity preserves the Weierstrass equation literally. The second
preserves it after \(x\mapsto x/t^2\), \(y\mapsto y/t^3\). Thus \(E\) is
\(S_3\)-equivariant over \(\mathbb Q(t)\).

## The three characters

Choose

\[
D_1=-t,
\qquad D_2=t-1,
\qquad D_3=-t(t-1)=D_1D_2.
\]

Their branch union is \(\{0,1,\infty\}\), so the associated biquadratic base
has genus zero. Under \(s\), \(D_1\leftrightarrow D_2\) and \(D_3\) is fixed.
Under \(i\), modulo squares,

\[
D_1\mapsto D_1/t^2,
\qquad D_2\mapsto D_3/t^2,
\qquad D_3\mapsto D_2/t^2.
\]

Therefore the three quadratic twists are pairwise isomorphic over
\(\mathbb Q(t)\).

## Theorem

Let

\[
r_0=\operatorname{rank}E(\mathbb Q(t)),
\qquad
r=\operatorname{rank}E^{(-t)}(\mathbb Q(t)).
\]

Character decomposition over the biquadratic extension gives

\[
\boxed{
\operatorname{rank}E\bigl(\mathbb Q(t)(\sqrt{-t},\sqrt{t-1})\bigr)
=r_0+3r.
}
\]

In particular, a split rank-eight member together with \(r\ge8\) gives a
rational genus-zero family of rank at least

\[
8+3\cdot8=32.
\]

The first forced packet rank above the current record is therefore 32 rather
than 30. One high-rank twist is enough; the other two channels are forced by
symmetry.

## Why this changes the search

Without symmetry, a rank-thirty packet needs three twist ranks whose sum is at
least 22. Their Noether--Lefschetz conditions appear heavily overdetermined.
The \(S_3\) action makes the three rank conditions identical, reducing the
construction to two exact tasks:

1. find a member whose geometric \(E_8\) lattice descends completely to
   \(\mathbb Q(t)\);
2. prove that the single representative twist by \(-t\) has rank at least
   eight over \(\mathbb Q(t)\).

The current Groebner and finite-field workstreams test these tasks. This
compression theorem does not claim that the two rank conditions have already
been met simultaneously.
