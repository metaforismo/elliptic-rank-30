# Genus-zero three-character quadratic rank packet

## The family

Let

\[
K=\mathbf Q(t),\qquad
E_t:\ y^2=f_t(x)
=x^3-\frac{t^2+1}{t+1}x+t.
\]

Put

\[
d_1=t,
\qquad
d_2=\frac{2t}{t+1},
\qquad
L=K(u,v),\quad u^2=d_1,\quad v^2=d_2.
\]

There are three immediate points over \(L\):

\[
P_u=(0,u),\qquad
P_v=(1,v),\qquad
P_{uv}=(-1,uv),
\]

because the following identities hold exactly in \(\mathbf Q(t)\):

\[
f_t(0)=d_1,
\qquad
f_t(1)=d_2,
\qquad
f_t(-1)=d_1d_2=\frac{2t^2}{t+1}.
\]

Equivalently, if the standard short-model quadratic twist by \(d\) is written

\[
E_t^{(d)}:\ Y^2=X^3+d^2a(t)X+d^3b(t),
\]

then the three character twists possess \(K\)-rational sections.  The
executable certificate verifies these twist equations directly rather than
assuming a twist convention.

## Theorem

The extension \(L/K\) is a connected biquadratic cover of degree four and
its function field is rational over \(\mathbf Q\).  Moreover,

\[
\boxed{
\operatorname{rank}E_t(L)
\ge
\operatorname{rank}E_t(K)+3.
}
\]

The result is unconditional.

## Proof

### 1. Squareclass independence and genus

At the places \(t=0\), \(t=-1\), and \(t=\infty\), the valuation-parity
vectors are

\[
d_1:(1,0,1),
\qquad
d_2:(1,1,0),
\qquad
d_1d_2:(0,1,1).
\]

The first two vectors have rank two over \(\mathbf F_2\).  Hence \(d_1,d_2\)
are independent squareclasses and \([L:K]=4\).  Their total branch support is
exactly \(\{0,-1,\infty\}\).  Every branch point has inertia group of order
two, so Riemann--Hurwitz gives

\[
2g(L)-2
=4(-2)+3\cdot 4\left(1-\frac12\right)
=-2.
\]

Thus \(g(L)=0\).

An explicit rational parameter makes this conclusion constructive.  Set
\(w=v/u\).  Then

\[
v^2+w^2=2.
\]

The conic is parameterized by

\[
v=\frac{r^2-2r-1}{r^2+1},
\qquad
w=\frac{1-2r-r^2}{r^2+1},
\qquad
u=\frac vw,
\qquad
t=u^2.
\]

Consequently \(L\cong\mathbf Q(r)\), and the three points become explicit
rational sections over \(\mathbf Q(r)\):

\[
(0,u),\qquad(1,v),\qquad(-1,uv).
\]

The certificate checks the polynomial identities
\(v^2+w^2=2\), \(t=u^2\), and \(2t/(t+1)=v^2\) by exact coefficient
arithmetic.

### 2. Each character section is non-torsion

Specialize the three twist sections at the smooth fiber \(t=-3\).  The base
curve is

\[
y^2=x^3+5x-3,
\qquad \Delta=-11888,
\]

and the three twist parameters are \(-3,3,-9\).  In integral short
Weierstrass models the specialized points are

\[
\begin{array}{c|c|c}
d&E_{-3}^{(d)}&Q_d\\
\hline
-3&y^2=x^3+45x+81&(0,9)\\
3&y^2=x^3+45x-81&(3,9)\\
-9&y^2=x^3+405x+2187&(9,81).
\end{array}
\]

Exact group-law arithmetic gives

\[
2Q_{-3}=\left(\frac{25}{4},-\frac{197}{8}\right),
\]

\[
3Q_3=\left(\frac{1479}{49},\frac{58185}{343}\right),
\]

and

\[
3Q_{-9}=\left(\frac{13077}{121},-\frac{1522395}{1331}\right).
\]

If one of the \(Q_d\) were torsion, the displayed nonzero multiple would also
be torsion.  Lutz--Nagell would then force integral coordinates on its
integral short Weierstrass model, contradicting the nonintegral
\(x\)-coordinate.  Therefore all three specialized points are non-torsion.
A torsion section specializes to a torsion point at every smooth fiber where
it is defined, so the three generic twist sections are non-torsion as well.

### 3. Character separation proves independence

Let \(G=\operatorname{Gal}(L/K)=\langle\sigma_u,\sigma_v\rangle\), where the
two generators change the signs of \(u\) and \(v\).  The sign characters of
the points are

\[
\begin{array}{c|cc}
&\sigma_u&\sigma_v\\
\hline
P_u&-1&+1\\
P_v&+1&-1\\
P_{uv}&-1&-1.
\end{array}
\]

After tensoring \(E_t(L)\) with \(\mathbf Q\), the commuting involutions split
the group into character eigenspaces.  The idempotent

\[
e_\chi=\frac14\sum_{g\in G}\chi(g)g
\]

projects onto the \(\chi\)-eigenspace.  Applying the three nontrivial
projectors to any proposed relation isolates the corresponding point.  Since
each point is non-torsion, each is nonzero in \(E_t(L)\otimes\mathbf Q\), and
the three points are independent.  The trivial-character projector also
shows that these directions are independent from \(E_t(K)\).  This proves the
rank inequality.

## Research significance and limitation

This theorem exhibits a structural mechanism that is easy to miss when two
quadratic multisections are counted separately: their product character can
carry a third independent section without increasing the genus of the
composite cover beyond zero.

For a high-generic-rank family, the correct audit is therefore not merely
"how many quadratic multisections were constructed?" but:

1. determine the complete squareclass span of the multisections;
2. inspect every nontrivial character twist, including product characters;
3. prove non-torsion of every surviving character section;
4. compute the total branch support and auxiliary genus;
5. rank constructions by the total character-packet contribution.

The displayed family itself certifies three generic directions; it is **not**
a rank-30 curve.  Its value is as a construction module to transplant into a
family already carrying many independent sections.

The executable proof is `research/packet3_certificate.py`.  Its deterministic
machine-readable output is `certificates/packet3_certificate.json`.
