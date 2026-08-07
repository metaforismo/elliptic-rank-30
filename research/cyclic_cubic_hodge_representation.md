# The C3 Hodge representation of a cubic genus-one pullback

**Author:** Francesco Giannicola  
**Truth status:** new intermediate theorem; no rank-30 curve is claimed.

## Setting

Let `S -> P1` be a rational elliptic surface with twelve `I1` fibres and split
Mordell--Weil lattice `E8`. Let

\[
f:C\to\mathbf P^1
\]

be a cyclic cubic cover of genus one, totally ramified over three smooth
fibres, and let

\[
S_C=S\times_{\mathbf P^1}C.
\]

Let `sigma` generate the deck group `C3`.

## Hodge numbers

The base-changed elliptic surface has

\[
\chi(\mathcal O_{S_C})=3,
\qquad g(C)=1.
\]

Therefore

\[
q=1,
\qquad p_g=\chi+g-1=3,
\]

and

\[
e=12\chi=36.
\]

Since

\[
e=2-4q+b_2,
\]

we obtain

\[
\boxed{b_2=38.}
\]

Moreover,

\[
\boxed{h^{1,1}=b_2-2p_g=32.}
\]

## Topological C3 representation

The fixed locus of `sigma` consists of the three smooth elliptic fibres above
the ramification points. Its Euler characteristic is zero. The topological
Lefschetz number is therefore

\[
L(\sigma)=0.
\]

On `H1`, the invariant part vanishes because `C/C3=P1`; the two eigenvalues
are `zeta3` and `zeta3^2`, so

\[
\operatorname{Tr}(\sigma|H^1)=-1.
\]

By Poincare duality the same trace occurs on `H3`. Hence

\[
0=L(\sigma)
=1-(-1)+\operatorname{Tr}(\sigma|H^2)-(-1)+1,
\]

and

\[
\operatorname{Tr}(\sigma|H^2)=-4.
\]

Write the complex eigenspace dimensions as

\[
\dim H^2_1=n_0,
\qquad
\dim H^2_{\zeta_3}=\dim H^2_{\zeta_3^2}=n.
\]

Then

\[
n_0+2n=38,
\qquad
n_0-n=-4.
\]

Thus

\[
\boxed{(n_0,n,n)=(10,14,14).}
\]

The invariant ten-dimensional space is exactly the pullback of `H2(S)`.

## Holomorphic two-forms

Let `M=O_P1(1)` define the cyclic cover, so `M^3=O(B)` for the degree-three
branch divisor. Then

\[
f_*O_C=O\oplus O(-1)\oplus O(-2)
\]

with the three deck characters. The canonical generator of the genus-one
cover contributes a nontrivial character. After the character shift, the
holomorphic two-forms have multiplicities

\[
\boxed{
H^{2,0}_1=0,
\qquad
(\dim H^{2,0}_{\zeta_3},\dim H^{2,0}_{\zeta_3^2})=(2,1)
}
\]

up to interchanging the two nontrivial characters.

Complex conjugation swaps the characters and the Hodge types. Therefore each
nontrivial `H2` eigenspace contains three non-`(1,1)` dimensions in total.
Consequently

\[
\boxed{
\dim H^{1,1}_1=10,
\qquad
\dim H^{1,1}_{\zeta_3}
=
\dim H^{1,1}_{\zeta_3^2}
=11.
}
\]

## Rank-thirty equivalence

Because every singular fibre remains irreducible, Shioda--Tate gives

\[
\rank\MW(S_C)=\rho(S_C)-2.
\]

The invariant Mordell--Weil lattice is the pulled-back `E8(3)` of rank eight.
The nontrivial character spaces can contribute at most eleven directions
each. Hence

\[
\boxed{
\rank\MW(S_C)\le8+11+11=30.
}
\]

The following are equivalent:

\[
\boxed{
\begin{aligned}
&\rank\MW(S_C)=30,\\
&\rho(S_C)=32=h^{1,1}(S_C),\\
&\dim\MW_{\zeta_3}=\dim\MW_{\zeta_3^2}=11,\\
&\text{every nontrivial-character }(1,1)\text{ class is algebraic.}
\end{aligned}
}
\]

Thus a rank-thirty cubic pullback is a Picard-maximal irregular elliptic
surface. Its trace-zero Mordell--Weil lattice is an Eisenstein lattice of
complex rank eleven.

## Consequence for search design

The target `11` is not an arbitrary count chosen to complement the base rank
eight. It is the full nontrivial-character `H11` capacity. The trisection
search is therefore trying to algebraize an entire eleven-dimensional
Eisenstein Hodge space.

This connects the constructive problem to:

* Picard-maximal surfaces;
* Noether--Lefschetz loci with an order-three automorphism;
* CM and Delsarte specializations;
* Eisenstein lattices and ternary trace codes;
* outer Galois points on marked cubic surfaces.

A successful packet cannot leave even one nontrivial `(1,1)` direction
transcendental. This explains why rank thirty should occur only at highly
special arithmetic points while still leaving a finite constructive route.
