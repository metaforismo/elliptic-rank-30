# Six character channels and the exact rank-30 capacity criterion

**Author:** Francesco Giannicola  
**Truth status:** new intermediate theorem; rational descent and explicit
rank-30 sections remain open.

## Generic marked family

Let

\[
q(t)=t^2-t+1,\qquad r(t)=t(t-1),
\]

choose a rational marking \(\mu\notin\{0,1\}\), and put

\[
a=\mu^2-\mu+1,
\qquad
b=\frac{a^3}{\mu(\mu-1)}.
\]

On the auxiliary coordinate

\[
v=t(t-1),
\]

the marked \(j=0\) surface is

\[
E_v:\quad y^2=x^3+P_\mu(v),
\qquad
P_\mu(v)=b^2v^2-a^3(v+1)^3.
\]

The cyclic genus-one base is obtained by adjoining

\[
u^3=cv,
\qquad
\sqrt{1+4v}=2t-1.
\]

Over an algebraic closure the constant \(c\) is a cube, so the deck group is
\(C_3\times C_2\).

## Rational splitting of the marked coefficient

Exact factorization gives

\[
\boxed{
P_\mu(v)=
-\frac{a^3}{\mu^2(\mu-1)^2}
\bigl(v-\mu(\mu-1)\bigr)
\bigl(\mu^2v+\mu-1\bigr)
\bigl((\mu-1)^2v-\mu\bigr).
}
\]

Thus its three roots are rational:

\[
r_1=\mu(\mu-1),
\qquad
r_2=\frac{1-\mu}{\mu^2},
\qquad
r_3=\frac{\mu}{(\mu-1)^2}.
\]

They satisfy

\[
r_1r_2r_3=-1,
\qquad
r_1r_2+r_1r_3+r_2r_3=3.
\]

The roots are distinct and avoid \(0\) and \(-1/4\) precisely away from

\[
\mu\in\left\{0,1,\frac12,2,-1\right\}.
\]

These five markings are the degenerate locus. The previously studied
\(\mu=2\) family lies on it, which explains its repeated factors and rank
collapse.

## Six-channel decomposition

For a \(j=0\) curve, a cubic twist by \(v^i\) multiplies \(a_6\) by
\(v^{2i}\), and a quadratic twist by \(D=1+4v\) multiplies \(a_6\) by
\(D^3\). After tensoring with \(\mathbf Q(\zeta_3)\), the Mordell--Weil
space over the full \(C_2\times C_3\) extension decomposes into the six
surfaces

\[
E_{i,\varepsilon}:\quad
 y^2=x^3+v^{2i}(1+4v)^{3\varepsilon}P_\mu(v),
\qquad
 i=0,1,2,
\quad \varepsilon=0,1.
\]

## Exact capacity vector

Assume \(\mu\) is outside the degenerate locus. The three roots of
\(P_\mu\), the point \(v=0\), the point \(v=-1/4\), and infinity are all
distinct. Kodaira classification and Shioda--Tate give the following exact
capacity table.

| \((i,\varepsilon)\) | \(\chi\) | additional reducible fibres | root rank | rank capacity |
|---|---:|---|---:|---:|
| \((0,0)\) | 1 | \(I_0^*\) at infinity | 4 | 4 |
| \((1,0)\) | 1 | \(IV\) at 0 | 2 | 6 |
| \((2,0)\) | 2 | \(IV^*\) at 0, \(II^*\) at infinity | 14 | 4 |
| \((0,1)\) | 1 | \(I_0^*\) at \(-1/4\) | 4 | 4 |
| \((1,1)\) | 2 | \(IV\), \(I_0^*\), \(IV^*\) | 12 | 6 |
| \((2,1)\) | 2 | \(IV^*\), \(I_0^*\), \(IV\) | 12 | 6 |

The capacity vector is therefore

\[
\boxed{(4,6,4,4,6,6)}
\]

and its sum is exactly

\[
\boxed{30}.
\]

The three \(\chi=1\) channels are rational elliptic surfaces, so their
geometric Picard number is automatically ten. Their geometric Mordell--Weil
ranks are therefore exactly

\[
4+6+4=14.
\]

The remaining three channels are K3 surfaces with rank capacities

\[
\boxed{4,6,6}.
\]

## Picard-maximality criterion

The generic cubic packet reaches geometric rank 30 if and only if all three
K3 channels are simultaneously singular K3 surfaces:

\[
\boxed{
\rho(E_{2,0})=\rho(E_{1,1})=\rho(E_{2,1})=20.
}
\]

Indeed the rational-surface channels already contribute 14, while the K3
channels can contribute at most 16. Equality in the total Hodge bound is
therefore equivalent to equality in every K3 Picard bound.

This replaces an unstructured search for 22 extra sections by three exact
Picard-maximality tests.

## Arithmetic requirements

Geometric rank 30 is still not a rank-30 curve over \(\mathbf Q\). A
constructive success additionally requires:

1. the six character lattices to descend through the rational cubic and
   quadratic twists;
2. a positive-rank genus-one base \(u^3=cv\);
3. explicit sections in all three K3 channels;
4. a rational specialization preserving 30 independent directions;
5. the complete height, reduction, saturation, and independent-CAS
   certificate.

The scale \(c\) does not affect geometric capacity, but it changes rational
descent and the Mordell--Weil rank of the genus-one base. Exact Sage
calculations already identify positive-rank scale classes such as
\(c=3,18,24,30,36,81\) in the tested range.

## Cheapest decisive experiment

For a rational marking \(\mu\), reduce the three K3 channels modulo several
good primes and compute their geometric function-field ranks. If any channel
has rank below its target \((4,6,6)\) at one good prime, the characteristic-zero
marking cannot be Picard maximal and is rejected.

Only markings that attain \((4,6,6)\) at every completed prime are promoted
to characteristic-zero section construction. The repository implements this
as the official-Magma six-channel sieve.
