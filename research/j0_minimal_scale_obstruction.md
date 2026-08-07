# Minimal-character scale obstruction on the marked j=0 cubic packet

**Author:** Francesco Giannicola  
**Truth status:** restricted-family theorem; no rank-30 curve is claimed.

## Setting

Let

\[
q(t)=t^2-t+1,\qquad r(t)=t(t-1)=q(t)-1,
\]

and take the marked parameter \(\mu=2\). Then

\[
a=q(2)=3,\qquad b=\frac{q(2)^3}{r(2)}=\frac{27}{2}.
\]

The rational elliptic surface is

\[
E:\quad y^2=x^3+b^2r(t)^2-a^3q(t)^3.
\]

For a rational nonzero scale \(c\), make the cyclic cubic base change

\[
C_c:\quad u^3=c\,t(t-1).
\]

Over \(K=\mathbf Q(w)\), \(w^2=-3\), the pullback becomes

\[
y^2=x^3+\frac{b^2}{c^2}u^6-a^3\left(\frac{u^3}{c}+1\right)^3.
\]

We classify every character-one section in the minimal polynomial box

\[
x=u(Au^3+B),\qquad y=Cu^6+Du^3+E,
\qquad A,B,C,D,E\in K.
\]

## Theorem

For the surface above, a minimal character-one eigensection exists over
\(K(C_c)\) only when

\[
[c]\in\{[1],[2]\}\subset \mathbf Q^*/\mathbf Q^{*3}.
\]

Both classes occur, and each carries a non-torsion CM eigenline over
\(K(C_c)\). However, the two corresponding genus-one bases \(C_1\) and
\(C_2\) both have exact Mordell--Weil rank zero over \(\mathbf Q\).
Consequently no positive-rank base in this \(\mu=2\) family can carry a
minimal character-one eigensection.

This does not exclude controlled-pole sections, another marking parameter,
another branch triple, or another elliptic surface.

## Complete coefficient reduction

Coefficient comparison gives

\[
\begin{aligned}
C^2&=A^3,\\
2CD&=3A^2B-\frac{27}{c^3},\\
D^2+2CE&=3AB^2+\frac{405}{4c^2},\\
2DE&=B^3-\frac{81}{c},\\
E^2&=-27.
\end{aligned}
\]

The branch \(B=0\) is impossible. Indeed the last and fourth equations give
\(D=3E/(2c)\), the second then gives \(C=E/(3c^2)\), and hence

\[
D^2+2CE=-\frac{315}{4c^2},
\]

contrary to the required \(405/(4c^2)\). Thus \(B\ne0\).

Choose \(E=3w\) and define

\[
\rho=\frac{cA}{B},\qquad
S=cB^3,\qquad g=\frac{c^2C}{w},\qquad d=\frac{cD}{w}.
\]

The five equations reduce to

\[
\begin{aligned}
-3g^2&=\rho^3S,\\
-2gd&=\rho^2S-9,\\
-3d^2-18g&=3\rho S+\frac{405}{4},\\
S&=81-18d.
\end{aligned}
\]

Exact lexicographic Groebner elimination yields

\[
\boxed{
(2\rho+1)^3
\left(10\rho^3-48\rho^2+12\rho-1\right)=0.
}
\]

The cubic factor is irreducible over \(\mathbf Q\) and has discriminant
\(-78732\). A root of an irreducible cubic cannot lie in the quadratic field
\(K\). Hence

\[
\rho=-\frac12.
\]

Substitution gives exactly two normalized solutions:

\[
(S,g,d)=\left(54,-\frac32,\frac32\right),
\qquad
(S,g,d)=\left(216,3,-\frac{15}{2}\right).
\]

## Kummer descent of the scale

If \(z\in K\) and \(z^3\in\mathbf Q\), then \(z^3\in\mathbf Q^3\). To see
this, write \(z=x+yw\). The coefficient of \(w\) in \(z^3\) is
\(3y(x^2-y^2)\), so either \(y=0\) or \(x=\pm y\); in every case \(z^3\)
is a rational cube. Thus

\[
K^3\cap\mathbf Q=\mathbf Q^3.
\]

Since \(S/c=B^3\in K^3\cap\mathbf Q\), the rational cube class of \(c\)
equals that of \(S\). Now

\[
54=2\cdot3^3,\qquad216=6^3,
\]

so the only possible classes are \([2]\) and \([1]\).

## Explicit eigenlines

For \(c=1\), one section is

\[
\begin{aligned}
x&=u(-3u^3+6),\\
y&=3wu^6-\frac{15}{2}wu^3+3w.
\end{aligned}
\]

For \(c=2\), one section is

\[
\begin{aligned}
x&=u\left(-\frac34u^3+3\right),\\
y&=-\frac38wu^6+\frac34wu^3+3w.
\end{aligned}
\]

Exact substitution verifies both identities.

At \(u=1\), finite reductions certify non-torsion. For the class \([1]\),
the reduced point has orders \(13\) at a prime above \(7\) and \(19\) at a
prime above \(13\). A torsion order would have to satisfy simultaneously

\[
n=13\cdot7^a=19\cdot13^b,
\]

which is impossible. For the class \([2]\), the corresponding orders are
\(13\) and \(21\), forcing

\[
n=13\cdot7^a=21\cdot13^b,
\]

again impossible.

The CM orbit \(P,[\zeta_3]P,[\zeta_3^2]P\) satisfies

\[
P+[\zeta_3]P+[\zeta_3^2]P=O.
\]

Since \(P\) is non-torsion, \(P\) and \([\zeta_3]P\) are independent over
\(\mathbf Z\), giving one Eisenstein line, or two ordinary lattice
directions, over \(K(C_c)\).

## Exact base-rank incompatibility

The base curve \(C_c\) is birational to

\[
Y^2=X^3+16c^4,
\qquad
X=4cu,\quad Y=4c^2(2t-1).
\]

SageMath 10.9 with `proof.all(True)` and `rank(proof=True)` gives

\[
\operatorname{rank}C_1(\mathbf Q)=0,
\qquad
\operatorname{rank}C_2(\mathbf Q)=0,
\]

and torsion of order three in both cases. Therefore the only minimal
eigensection-compatible cube classes do not provide infinitely many rational
specialization parameters.

## What the failure reveals

The earlier computations showed an apparent tradeoff:

- \(c=27\), which is cube-equivalent to \(c=1\), carries explicit minimal
  eigensections but its base has rank zero;
- \(c=36\) has base rank one but the minimal eigensection scheme is empty.

The theorem proves this is not an accident of those two representatives. In
the entire \(\mu=2\) minimal polynomial box, geometric visibility and rational
base rank are disjoint.

## Decisive change of setting

Do not enlarge the same coefficient box. The next search must do at least one
of the following:

1. allow controlled poles in the character-one sections on a positive-rank
   cube class such as \(c=3\) or \(c=36\);
2. move along the CM-compatible marking conic
   \(q(\mu)=3s^2\) and repeat the exact Kummer classification;
3. vary the rational branch triple rather than fixing \(0,1,\infty\);
4. return to the full split-\(E_8\) trace-code incidence search.

The cheapest next experiment is the one-pole character-one ansatz on the
rank-one base \(C_3\).
