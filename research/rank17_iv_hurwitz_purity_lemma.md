# Hurwitz-purity lemma for the normalized additive-\(IV\) locus

## Statement

Let \(U_{IV}\) be the characteristic-zero open locus of normalized pairs

\[
c_4=(t-1)^2A(t),\qquad c_6=(t-1)^2B(t),
\]

where \(A,B\) are monic of degrees \(6,10\), and

\[
(t-1)^2A(t)^3-B(t)^2=\kappa\,t^4R_4(t)
\]

with \(\kappa\ne0\) and \(R_4\) squarefree of degree \(4\). Impose the usual
coprimality and noncollision conditions:

- \(A,B,R_4\) have the required simple roots;
- \(t=0,1,\infty\) do not collide with the other marked points;
- the extra critical point is simple;
- its critical value is different from \(0,1,\infty\).

Then every complex point of \(U_{IV}\) has local dimension one. In
particular, \(U_{IV}\) has no isolated geometric points and no
zero-dimensional irreducible component.

This statement concerns the open nondegenerate surface locus. It says
nothing about boundary points where a fibre worsens, the discriminant
vanishes identically, or two branch values collide.

## 1. The associated degree-20 cover

Define

\[
\phi(t)=
\frac{(t-1)^2A(t)^3}
     {(t-1)^2A(t)^3-B(t)^2}.
\]

The three fixed fibres have passports

\[
\begin{aligned}
\phi^{-1}(0)&:\quad (3^6,2),\\
\phi^{-1}(1)&:\quad (2^{10}),\\
\phi^{-1}(\infty)&:\quad (12,4,1^4).
\end{aligned}
\]

Indeed,

\[
\phi-1=
\frac{B(t)^2}
     {(t-1)^2A(t)^3-B(t)^2},
\]

and the pole of order \(12\) at \(t=\infty\) is the difference between the
degree \(20\) numerator and the degree \(8\) affine denominator.

The ramification contributions at the three fixed branch values are

\[
13,\qquad 10,\qquad 14.
\]

Their sum is \(37\). Since

\[
2\deg(\phi)-2=38,
\]

there is exactly one further simple ramification point. Its branch cycle is

\[
(2,1^{18}).
\]

Thus the ordered passport is

\[
\boxed{
(3^6,2),\quad
(2^{10}),\quad
(12,4,1^4),\quad
(2,1^{18}).
}
\]

The fourth branch value will be denoted by \(\lambda\).

## 2. The Hurwitz branch map is locally unramified

For covers of the Riemann sphere with a fixed Nielsen class, the Hurwitz
space maps to the configuration space of branch points as an unramified
covering. Fried and Völklein construct precisely this topology and state
that the branch maps are unramified coverings; the Hurwitz spaces thereby
inherit complex-manifold structures.

Primary reference:

- Michael D. Fried and Helmut Völklein,
  *The inverse Galois problem and rational points on moduli spaces*,
  **Math. Ann. 290** (1991), 771–800, §1.2, especially pp. 774–775.
  DOI: `10.1007/BF01459271`.
- Their construction refers to William Fulton,
  *Hurwitz schemes and irreducibility of moduli of algebraic curves*,
  **Ann. of Math. 90** (1969), 542–575, §1.3.
  DOI: `10.2307/1970748`.

After ordering the four branch points and fixing the first three at

\[
0,\quad1,\quad\infty,
\]

the branch configuration space is

\[
\mathbf P^1_\lambda\setminus\{0,1,\infty\}.
\]

It is a smooth complex curve. Therefore every Hurwitz-space point with
this passport has a one-dimensional analytic neighbourhood.

## 3. Why our coefficient normalization is a genuine local chart

The source coordinate is rigidified by three uniquely characterized
ramification points:

- \(t=1\) is the unique point of index \(2\) above branch value \(0\);
- \(t=0\) is the unique point of index \(4\) above branch value \(\infty\);
- \(t=\infty\) is the unique point of index \(12\) above branch value
  \(\infty\).

Any automorphism of the source preserving the normalized cover fixes these
three points, hence is the identity.

The target coordinate is rigidified by the ordered values
\(0,1,\infty\). The monicity of \(A\) and \(B\) fixes the remaining cube-
and square-root scalars. Consequently, near every nondegenerate point,
the normalized coefficient description and the corresponding ordered
Hurwitz space determine one another uniquely.

Thus the normalized coefficient locus is locally analytically isomorphic
to an open subset of the reduced Hurwitz space over the \(\lambda\)-line.

## 4. Consequence for isolated algebraic solutions

Analytic and algebraic local dimensions agree for a finite-type complex
scheme. Hence every point of \(U_{IV}(\mathbf C)\) lies on a
one-dimensional algebraic component.

In particular, a rational point of the nondegenerate open locus cannot be
an isolated solution of the residual polynomial system. After embedding

\[
\mathbf Q\hookrightarrow\mathbf C,
\]

it must lie on one of the positive-dimensional components classified by
the exact resultant computation.

## Scope

This lemma is used only after verifying the open conditions that give the
four distinct branch values and the passport above. It cannot be applied
to:

- \(\kappa=0\);
- repeated roots of the residual quartic;
- collisions \(\lambda\in\{0,1,\infty\}\);
- worsening of the fibres at \(0,1,\infty\);
- points at infinity of a coefficient chart that do not represent a
  normalized elliptic K3 surface.
