# Torsion forces too much reducible-fibre cost for a rank-30 smooth packet

Let `S -> P^1` be a rational elliptic surface over a field of characteristic
zero.  Write `r0` for its geometric Mordell--Weil rank.  For a connected split
`V4` base change ramified at three smooth fibres, let `r1,r2,r3` be the three
quadratic-twist ranks.  The smooth-branch packet bound is

\[
r_0+r_1+r_2+r_3\le 4r_0+6.
\]

The following elementary fibre-cost bounds rule out the most tempting torsion
models.

## Rational two-torsion

After choosing the two-torsion section as `(0,0)`, a minimal global model has

\[
y^2=x^3+A(t)x^2+B(t)x,
\qquad A\in H^0(\mathbb P^1,\mathcal O(2)),
\quad B\in H^0(\mathbb P^1,\mathcal O(4)).
\]

Its discriminant is

\[
\Delta=16B^2(A^2-4B).
\]

If `ord_v(B)=m>0`, Tate's algorithm gives a fibre-root contribution at least
`2m-1`: it is `I_{2m}` when `A` is a unit, and the additive collision cases
have at least the same root-lattice rank.  Since the zero divisor of `B` has
degree four,

\[
\sum_v(m_v-1)\ge
\sum_{B(v)=0}(2\operatorname{ord}_v(B)-1)
=8-\#\operatorname{Supp}(B)\ge4.
\]

Shioda--Tate on a rational elliptic surface therefore gives

\[
r_0\le 8-4=4.
\]

Consequently every smooth three-branch packet built from such a surface obeys

\[
\boxed{r_0+r_1+r_2+r_3\le22.}
\]

It cannot produce generic rank 30.

## Rational three-torsion

In Tate normal form,

\[
y^2+a_1xy+a_3y=x^3,
\qquad a_1\in H^0(\mathbb P^1,\mathcal O(1)),
\quad a_3\in H^0(\mathbb P^1,\mathcal O(3)),
\]

and

\[
\Delta=a_3^3(a_1^3-27a_3).
\]

At a zero of `a3` of multiplicity `m`, the root contribution is at least
`3m-1`; the generic case is `I_{3m}` and the additive collision cases meet the
same lower bound.  Since `div_0(a3)` has degree three,

\[
\sum_v(m_v-1)\ge 9-\#\operatorname{Supp}(a_3)\ge6,
\]

so

\[
r_0\le2
\quad\text{and}\quad
\boxed{r_0+r_1+r_2+r_3\le14.}
\]

## Search consequence

A certificate-first rank-30 packet search should not prioritize rational
2-torsion or 3-torsion merely because their descent maps are easy to write.
The reducible fibres consume more rank capacity than the simplified local
conditions can repay.  The high-capacity rational base must have trivial small
rational torsion and preferably twelve irreducible `I1` fibres.

This is a restricted obstruction to the smooth three-branch `V4` mechanism;
it is not a universal upper bound for elliptic curves over `Q`.
