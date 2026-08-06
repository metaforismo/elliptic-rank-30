# Genus-zero three-channel quadratic rank packet

## The family

Let

\[
K=\mathbb Q(t),\qquad
E_t:\ y^2=x^3-\frac{t^2+1}{t+1}x+t.
\]

Put

\[
d_1=t,
\qquad d_2=\frac{2t}{t+1}.
\]

The three nontrivial character twists in the biquadratic extension

\[
L=K(\sqrt{d_1},\sqrt{d_2})
\]

contain the following sections:

\[
\begin{array}{c|c|c}
\text{twist}&x&y_{\rm twist}\\
\hline
d_1y^2=f_t(x)&0&1\\
d_2y^2=f_t(x)&1&1\\
d_1d_2y^2=f_t(x)&-1&1,
\end{array}
\]

where

\[
f_t(x)=x^3-\frac{t^2+1}{t+1}x+t.
\]

Indeed,

\[
f_t(0)=t,
\qquad
f_t(1)=\frac{2t}{t+1},
\qquad
f_t(-1)=\frac{2t^2}{t+1}=d_1d_2.
\]

## Theorem

The squareclasses of \(d_1\) and \(d_2\) are geometrically independent. Their
total geometric branch support on \(\mathbb P^1_t\) is

\[
\{0,-1,\infty\}.
\]

Consequently, the connected biquadratic cover attached to
\(\langle d_1,d_2\rangle\) has genus zero. The three displayed sections are
non-torsion and occupy the three distinct nontrivial characters of
\(\operatorname{Gal}(L/K)\). Therefore

\[
\boxed{
\operatorname{rank}E_t(L)
\ge
\operatorname{rank}E_t(K)+3.
}
\]

### Proof

At the places \(0,-1,\infty\), the parity vectors of the two squareclasses are

\[
d_1:(1,0,1),
\qquad d_2:(1,1,0).
\]

They are linearly independent over \(\mathbb F_2\), and their sum is
\((0,1,1)\). Thus the geometric character dimension is \(2\) and the union of
branch points has size \(3\). Riemann--Hurwitz gives

\[
g=1+2^{2-2}(3-4)=0.
\]

To prove that the sections are non-torsion, specialize at \(t=-3\). The base
curve is

\[
y^2=x^3+5x-3,
\]

and the three twist parameters are \(-3,3,-9\). In standard integral twist
models the relevant points are

\[
\begin{array}{c|c|c}
d& E^{(d)} & P_d\\
\hline
-3&y^2=x^3+45x+81&(0,9)\\
3&y^2=x^3+45x-81&(3,9)\\
-9&y^2=x^3+405x+2187&(9,81).
\end{array}
\]

Exact addition gives

\[
x(2P_{-3})=\frac{25}{4},
\qquad
x(3P_3)=\frac{1479}{49},
\qquad
x(3P_{-9})=\frac{13077}{121}.
\]

If any \(P_d\) were torsion, all of its multiples would be torsion; by the
Lutz--Nagell theorem their coordinates on these integral short Weierstrass
models would be integral. The displayed nonintegral coordinates prove that all
three specialized points, and hence all three generic twist sections, are
non-torsion.

The Néron--Tate pairing is Galois invariant. Sections belonging to distinct
quadratic characters are mutually orthogonal, hence independent. This proves
the rank inequality.

## Explicit rational parametrisation

Write

\[
u^2=t,
\qquad
v^2=\frac{2t}{1+t},
\qquad
w=\frac vu.
\]

Then \(v^2+w^2=2\). The conic has the parametrisation

\[
v=\frac{r^2-2r-1}{r^2+1},
\qquad
w=\frac{1-2r-r^2}{r^2+1},
\qquad
u=\frac vw,
\qquad
t=u^2.
\]

Thus the biquadratic base is rational over \(\mathbb Q\), and the three points
become rational sections over \(\mathbb Q(r)\):

\[
(0,u),
\qquad
(1,v),
\qquad
(-1,uv).
\]

## Why this matters for rank 30

This gives an exact counterexample to the implicit accounting rule

\[
\text{two quadratic multisections}
\Longrightarrow
\text{at most two added directions}.
\]

The product character can carry a third independent section without any genus
penalty beyond the original composite cover. Therefore, for a rank-17 K3
fibration, a genus-zero pair with a positive product-twist channel forces rank
at least \(20\), not merely \(19\).

The decisive computational audit for every pair \((d_1,d_2)\) is now:

1. verify the two known twist sections;
2. compute the complete Mordell--Weil group of the hidden product twist
   \(E^{d_1d_2}\);
3. compute the union of branch supports and the auxiliary genus;
4. promote the pair by total packet rank, not by the number of multisections
   originally used to construct it.

The executable certificate is `research/packet3_certificate.py`; its
machine-readable output is `certificates/packet3_certificate.json`.
