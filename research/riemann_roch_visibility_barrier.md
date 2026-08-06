# The Riemann--Roch visibility barrier

Let `E/k` be an elliptic curve with identity `O`, and let `f` be a nonzero
rational function whose only pole is `O`, of exact order `N`.  If the zero
divisor consists of `N` rational points `P_1,...,P_N`, counted with
multiplicity, then

\[
(f)=P_1+\cdots+P_N-NO.
\]

Under `Pic^0(E) ~= E`, a principal divisor represents zero.  Therefore

\[
\boxed{P_1+\cdots+P_N=O.}
\]

The visible points from one function always have at least one exact integral
relation.

## Consequence for the near-square construction

On

\[
E:\ y^2=R(x),\qquad \deg R=3,
\]

let `Q(x)` have degree `n>=2`.  The function `f=y-Q(x)` has pole order `2n`
at `O`.  If

\[
Q(x)^2-R(x)=c\prod_{i=1}^{2n}(x-r_i)
\]

has `2n` distinct rational roots, then the points

\[
P_i=(r_i,Q(r_i))
\]

satisfy

\[
\sum_{i=1}^{2n}P_i=O.
\]

Thus a degree-30 split near-square can produce many rational points but can
never certify 30 independent points: its 30 visible points span rank at most
29.  This is a structural obstruction, not a failed search bound.

## The smallest one-function target that can still reach rank 30

Every function on a short Weierstrass curve can be written

\[
f=A(x)+yB(x).
\]

Its pole order is

\[
N=\max(2\deg A,\,2\deg B+3).
\]

The smallest `N` for which one principal-divisor relation is compatible with
rank 30 is `N=31`.  Take

\[
\deg A\le15,\qquad \deg B=14.
\]

Then

\[
\operatorname{Norm}_{k(E)/k(x)}(A+yB)
=A(x)^2-R(x)B(x)^2
\]

has degree 31.  A certificate-first construction target is

\[
\boxed{
A(x)^2-R(x)B(x)^2
=c\prod_{i=1}^{31}(x-r_i),\qquad r_i\in\mathbb Q.
}
\]

The corresponding 31 rational points have one forced relation, leaving room
for rank 30.  Independence still has to be proved; the norm identity alone is
not a rank certificate.

## Multi-function version

If rational points are assembled as zeros of several functions, every
independent principal divisor contributes a relation.  Candidate ranking must
therefore track the rank of the divisor-relation matrix, not merely the number
of visible points.  This explains why point counts and subgroup rank diverge
precisely at the frontier.
