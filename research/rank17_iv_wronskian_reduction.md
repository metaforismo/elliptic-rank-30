# Wronskian reduction of the additive-IV surface locus

## Status

This note proves an exact algebraic reduction of the normalized surface
problem. It does **not** construct a rational surface, a section of height
`79/12`, or a rank-30 elliptic curve.

## The almost-Belyi form

Write

\[
c_4(t)=(t-1)^2A(t),\qquad c_6(t)=(t-1)^2B(t),
\]

where `A` and `B` are monic of degrees `6` and `10`. Then

\[
\Delta_0(t)=c_4(t)^3-c_6(t)^2
=(t-1)^4\left((t-1)^2A(t)^3-B(t)^2\right).
\]

Put

\[
F=(t-1)^2A^3,\qquad G=B^2,\qquad H=F-G.
\]

The desired fibre orders are equivalent, on the relevant open subset, to

\[
H=t^4R_4(t),\qquad \deg R_4=4.
\]

The associated rational function

\[
\phi(t)=\frac{F(t)}{H(t)}
=\frac{(t-1)^2A(t)^3}{t^4R_4(t)}
\]

has the three distinguished ramification partitions

\[
3^6\,2,\qquad 2^{10},\qquad 12\,4\,1^4.
\]

Their ramification deficit is one. Thus the generic map has one further simple
ramification point: this is an almost-Belyi, one-dimensional Hurwitz problem.

## Exact Wronskian identity

Define

\[
C(t)=
2A(t)B(t)
+3(t-1)A'(t)B(t)
-2(t-1)A(t)B'(t).
\]

Direct differentiation gives

\[
\boxed{
F'G-FG'=(t-1)A^2B\,C.
}
\]

If `H=t^4R_4`, then also

\[
F'G-FG'=H'G-HG'.
\]

The right-hand side has the known factors from the pole of order four at zero
and the unique extra simple ramification point. Therefore

\[
\boxed{
C(t)=-12k\,t^3(t-e),
}
\]

where `k` is the leading coefficient of `R_4` and `e` is the extra
ramification point.

Equivalently, among the coefficients of `C`, only degrees `3` and `4` may be
nonzero.

## Converse

Assume characteristic zero, and

The finite-field replays below verify `H=t^4R_4` directly; they do not invoke this characteristic-zero degree argument.

\[
A(0)=p_0^2,\qquad B(0)=p_0^3,\qquad p_0\ne0.
\]

Suppose

\[
C(t)=c_3t^3+c_4t^4,\qquad c_3c_4B(1)\ne0.
\]

Because `F(0)=G(0)`, we have `H(0)=0`. If

\[
\operatorname{ord}_0(H)=m<4,
\]

then, since `G(0)=p_0^6` and `m=1,2,3` is nonzero in the field,

\[
\operatorname{ord}_0(H'G-HG')=m-1<3.
\]

This contradicts the factor `t^3` in the Wronskian. Since `c_3\ne0`, the order
is exactly four.

Next let `m=deg H`. The monic degree-20 terms of `F` and `G` cancel, and the
leading coefficient of

\[
H'G-HG'
\]

is

\[
(m-20)\operatorname{lc}(H)t^{m+19}.
\]

The factored Wronskian has degree `23+4=27`, because `c_4\ne0`. Hence

\[
m+19=27,
\]

so

\[
\boxed{\deg H=8.}
\]

Consequently

\[
\boxed{H=t^4R_4(t),\qquad \deg R_4=4.}
\]

Moreover,

\[
c_4=-12\operatorname{lc}(R_4),
\qquad
 e=-\frac{c_3}{c_4}.
\]

The condition `B(1) != 0` gives exact type `IV` at `t=1`. Nonzero `c_3` and
`c_4` give exact orders `I4` and `I12` at zero and infinity. Squarefreeness of
`R_4` and coprimality with `A` remain separate open conditions.

## Sparse equation system

Write

\[
A=p_0^2+a_1t+a_2t^2+a_3t^3+a_4t^4+a_5t^5+t^6,
\]

\[
B=p_0^3+b_1t+b_2t^2+\cdots+b_9t^9+t^{10}.
\]

The surface curve is cut out by the fourteen sparse equations

\[
[t^i]C=0,
\qquad
 i\in\{0,1,2,5,6,\ldots,15\}.
\]

Only `[t^3]C` and `[t^4]C` remain. The equations are linear in the `b_i` as a
block and mostly quadratic overall, replacing the dense degree-9 through
degree-11 square-root recurrence equations.

The first and last coefficients are already triangular. For example,

\[
[t^0]C=p_0^2(-3a_1p_0+2b_1+2p_0^3),
\]

\[
[t^{15}]C=-3a_5+2b_9+2.
\]

Thus one may eliminate `b_1,b_2,b_3` from the low end and
`b_9,b_8,\ldots,b_4` from the high end, leaving five equations in

\[
p_0,a_1,a_2,a_3,a_4,a_5.
\]

## Birational relation with the previous coordinates

Exact division of the original `c4` parametrization gives

\[
\begin{aligned}
a_1&=2p_0^2+2p_0p_1,\\
a_2&=3p_0^2+4p_0p_1+2p_0p_2+p_1^2,\\
a_3&=4p_0^2+6p_0p_1+4p_0p_2+2p_0p_3
     +2p_1^2+2p_1p_2,\\
a_4&=r+2s+3,\\
a_5&=s+2.
\end{aligned}
\]

Since `p0 != 0`, this change is birational. In particular,

\[
\boxed{s=a_5-2,\qquad r=a_4-2a_5+1.}
\]

Therefore an eliminant in `(a4,a5)` is immediately an eliminant in `(r,s)`.

## The exceptional F13 point

For the unique split surface found over `F13`, the Wronskian gives

\[
e=8.
\]

Its residual quartic is

\[
R_4(t)=-3(t+5)^2(t^2+2t+6).
\]

Thus `e=8=-5` is exactly the repeated residual pole. The fourth simple branch
point has collided with the pole fibre, explaining geometrically why this
point belongs to the broader split closure but not to the open configuration

\[
I_{12}+I_4+IV+4I_1.
\]
