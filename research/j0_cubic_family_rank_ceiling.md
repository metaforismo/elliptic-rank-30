# Hidden quadratic-twist decomposition and a rank-12 ceiling

**Author:** Francesco Giannicola  
**Truth status:** restricted-family theorem; special fibres may still jump.

## Statement

Let

\[
q(t)=t^2-t+1,\qquad r(t)=t(t-1),
\]

and take the marked parameter \(\mu=2\), so

\[
a=3,\qquad b=\frac{27}{2}.
\]

For every \(c\in\mathbf Q^*\), consider

\[
E_t:\quad y^2=x^3+b^2r(t)^2-a^3q(t)^3
\]

over the cyclic genus-one base

\[
C_c:\quad u^3=c\,t(t-1).
\]

Then

\[
\boxed{
\operatorname{rank}E(\overline{\mathbf Q}(C_c))\le12.
}
\]

Thus no member of this one-parameter cubic-base family can force generic rank
30. Any rank-30 specialization in it would require an exceptional jump of at
least 18.

## The hidden quadratic extension

Put

\[
v=t(t-1)=\frac{u^3}{c}.
\]

The elliptic-curve coefficient depends only on \(v\):

\[
\begin{aligned}
b^2v^2-a^3(v+1)^3
&=\frac{729}{4}v^2-27(v+1)^3\\
&=-\frac{27}{4}(v-2)^2(4v+1).
\end{aligned}
\]

Moreover \(t\) is quadratic over \(\mathbf Q(u)\), since

\[
t^2-t-\frac{u^3}{c}=0.
\]

Its discriminant is

\[
D_c(u)=1+\frac{4u^3}{c}=1+4v.
\]

Hence

\[
\mathbf Q(C_c)=\mathbf Q(u,\sqrt{D_c(u)}).
\]

For any elliptic curve over a field of characteristic zero, a quadratic
extension gives the rational-rank decomposition

\[
E(K(\sqrt D))\otimes\mathbf Q
\cong
E(K)\otimes\mathbf Q
\oplus
E^{D}(K)\otimes\mathbf Q.
\]

Therefore the geometric rank on \(C_c\) is the sum of the ranks of an
invariant surface over the \(u\)-line and its quadratic twist.

## Invariant surface

Over \(\overline{\mathbf Q}(u)\), the invariant curve is

\[
y^2=x^3-\frac{27}{4}
\left(\frac{u^3}{c}-2\right)^2
\left(\frac{4u^3}{c}+1\right).
\]

The finite fibres are:

- three fibres of type \(IV\) at \(u^3=2c\), contributing root rank
  \(3\cdot2=6\);
- three fibres of type \(II\) at \(4u^3=-c\), contributing no root lattice.

The coefficient has degree nine. In the minimal model of arithmetic genus
\(\chi=2\), the fibre at infinity has \(\operatorname{ord}(a_6)=3\), hence
type \(I_0^*\), contributing root rank four.

Thus

\[
R_{\mathrm{inv}}=6+4=10.
\]

The surface is a K3 surface, so

\[
h^{1,1}=20,
\qquad \rho\le20.
\]

Shioda--Tate gives

\[
\operatorname{rank}E_0(\overline{\mathbf Q}(u))
\le20-2-10=8.
\]

The Euler check is exact:

\[
3e(IV)+3e(II)+e(I_0^*)
=3\cdot4+3\cdot2+6=24=12\chi.
\]

## Quadratic-twist surface

For a short \(j=0\) model, the quadratic twist multiplies \(a_6\) by
\(D_c^3\). Therefore

\[
\begin{aligned}
a_6^{\mathrm{tw}}
&=D_c^3\left(-\frac{27}{4}(v-2)^2(4v+1)\right)\\
&=-\frac{27}{4}(v-2)^2(4v+1)^4.
\end{aligned}
\]

After substituting \(v=u^3/c\), the coefficient has degree 18, so the twist
surface has \(\chi=3\) and

\[
h^{1,1}=30.
\]

Its singular fibres are:

- three type-\(IV\) fibres at \(u^3=2c\), root rank six in total;
- three type-\(IV^*\) fibres at \(4u^3=-c\), root rank 18 in total;
- a smooth fibre at infinity.

Thus

\[
R_{\mathrm{tw}}=24
\]

and Shioda--Tate gives

\[
\operatorname{rank}E_0^{D_c}(\overline{\mathbf Q}(u))
\le30-2-24=4.
\]

Again the Euler check closes exactly:

\[
3e(IV)+3e(IV^*)=3\cdot4+3\cdot8=36=12\chi.
\]

Adding the two eigenspaces gives

\[
\boxed{
\operatorname{rank}E(\overline{\mathbf Q}(C_c))
\le8+4=12.
}
\]

## What the failure reveals

The degree-three genus-one Hodge bound of 30 is only a capacity bound. It can
be destroyed by reducible-fibre cost. In this symmetric marked family, the
coefficient factorization forces root-lattice rank

\[
10+24=34
\]

across the two quadratic eigenspaces, leaving room for only twelve
Mordell--Weil directions.

The earlier scale obstruction was therefore a symptom, not the whole disease.
Even controlled-pole sections cannot turn this particular family into a
generic rank-30 construction.

## Decisive change of setting

The primary cubic-cover search must move away from the special marking
\(\mu=2\), whose Legendre invariant is \(j=1728\) and whose coefficient
acquires repeated factors. The next family should satisfy all of the
following:

1. a split rank-eight rational elliptic surface before base change;
2. a positive-rank cyclic genus-one cubic base;
3. ramification over smooth fibres;
4. square-free or minimally reducible eigenspace coefficients after the
   hidden quadratic decomposition;
5. enough remaining Shioda--Tate capacity for an Eisenstein rank-eleven
   trace-zero lattice.

The cheapest next experiment is to run the same exact fibre-cost computation
symbolically along the CM-compatible conic \(q(\mu)=3s^2\), reject every
parameter whose two eigenspaces have combined rank ceiling below 30, and only
then construct trisection linear systems.
