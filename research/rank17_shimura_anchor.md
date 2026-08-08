# Exact Shimura anchor for the generic-rank-17 K3 surface

## Primary-source identification

Elkies's 2007 lectures identify the maximal-rank elliptic K3 construction with
a rational non-CM point on the genus-two Shimura curve

\[
X(6,79)/\langle w_{6\cdot79}\rangle.
\]

The same identification is repeated in the introduction of his 2008 paper on
Shimura-curve computations via K3 surfaces.

Primary sources:

- Noam D. Elkies, *Three lectures on elliptic surfaces and curves of high
  rank*, arXiv:0709.2908:
  <https://arxiv.org/abs/0709.2908>
- Noam D. Elkies, *Shimura curve computations via K3 surfaces of
  Néron-Severi rank at least 19*, arXiv:0802.1301:
  <https://arxiv.org/abs/0802.1301>

The published genus-two equation is

\[
C:\quad
u^2=16t^6-19t^4+88t^2-48.
\]

The source lists:

- the rational points at infinity fixed by the bielliptic involution;
- four affine points with \(|t|=2\), \(|u|=32\);
- four affine points with
  \(|t|=14/13\), \(|u|=2^6\cdot251/13^3\).

The latter ordinate is exactly

\[
\frac{2^6\cdot251}{13^3}
=
\frac{16064}{2197}.
\]

The source further states that the \(|t|=14/13\) orbit is non-CM and that it
yields an elliptic K3 surface over \(\mathbf Q(s)\) with Mordell-Weil rank 17.
Those moduli assertions are not yet independently reconstructed here.

## Exact facts certified in this repository

The executable certificate

```bash
python3 research/rank17_shimura_anchor.py \
  --compare certificates/rank17_shimura_anchor.json
```

uses only the Python standard library and proves the following.

### 1. The hyperelliptic curve is nonsingular of genus two

Let

\[
F(t)=16t^6-19t^4+88t^2-48.
\]

The exact resultant and discriminant are

\[
\operatorname{Res}(F,F')
=-960467703986795839488,
\]

\[
\operatorname{disc}(F)
=60029231499174739968\ne0.
\]

Thus \(F\) is squarefree. Since \(\deg F=6\), the smooth projective model of
\(u^2=F(t)\) has genus two.

The leading coefficient is \(16=4^2\), so the even-degree model has two
rational points at infinity. In weighted-projective coordinates
\((T:Z:U)\) of weights \((1,1,3)\), they are

\[
(1:0:4),\qquad(1:0:-4).
\]

### 2. The published affine rational points are exact

The certificate verifies

\[
F(2)=32^2
\]

and

\[
F\!\left(\frac{14}{13}\right)
=
\left(\frac{16064}{2197}\right)^2.
\]

Because \(F\) is even, this gives the eight distinct affine rational points

\[
(\pm2,\pm32)
\]

and

\[
\left(
\pm\frac{14}{13},
\pm\frac{16064}{2197}
\right),
\]

with all sign choices.

### 3. The displayed involution is genuinely bielliptic

The involution

\[
\iota(t,u)=(-t,-u)
\]

preserves \(C\). Its invariant functions

\[
x=t^2,\qquad y=tu
\]

satisfy

\[
y^2=16x^4-19x^3+88x^2-48x.
\]

The quartic discriminant is

\[
-80518053888\ne0,
\]

so the quotient is a nonsingular genus-one curve. On the open set
\(x\ne0\), the reciprocal change

\[
X=\frac1x=\frac1{t^2},
\qquad
Y=\frac{y}{x^2}=\frac{u}{t^3}
\]

gives the cubic model

\[
Y^2=-48X^3+88X^2-19X+16,
\]

whose cubic discriminant is

\[
-34947072\ne0.
\]

For example, the two positive representatives map to

\[
(t,u)=(2,32)
\longmapsto
(X,Y)=\left(\frac14,4\right)
\]

and

\[
\left(\frac{14}{13},\frac{16064}{2197}\right)
\longmapsto
\left(\frac{169}{196},\frac{2008}{343}\right).
\]

All orbit images are recorded exactly in the machine-readable certificate.

## What this anchor does not yet prove

The executable certificate does **not** prove any of the following source
claims:

1. the \(|t|=14/13\) orbit is non-CM;
2. its associated K3 surface has Mordell-Weil rank 17 over \(\mathbf Q(s)\);
3. that surface, in the relevant elliptic model, specializes to the published
   rank-29 curve.

The missing bridge is therefore precise:

\[
\boxed{
\text{Shimura point}
\longrightarrow
\text{QM/K3 moduli data}
\longrightarrow
\text{elliptic K3 equation}
\longrightarrow
17\text{ explicit sections}
\longrightarrow
E_{29}.
}
\]

## Construction route described by Elkies

The 2007 lectures explain the computational architecture used to cross this
bridge in the original work.

1. **Lattice design.** Construct the positive-definite essential lattice as a
   suitable slice of a Niemeier lattice, using roots and glue to control
   reducible fibers and torsion.
2. **Finite-field seed.** Solve the polynomial identities for the elliptic
   surface and sections modulo a suitable small prime.
3. **p-adic lifting.** Lift the finite-field solution by a multivariable
   Newton process, with finite differences approximating the Jacobian.
4. **Rational reconstruction.** Recognize the high-precision p-adic
   coefficients as rationals using lattice reduction, then verify every
   identity by exact substitution.
5. **Neighbor transformations.** Move between elliptic models of the same K3
   surface through chains of 2-neighbors and occasional 3-neighbors until the
   desired rootless essential lattice and manageable coefficients appear.
6. **Moduli deformation.** Starting from a known Picard-rank-20 point, deform
   p-adically along the Picard-rank-19 family, recover algebraic coefficient
   relations, verify them symbolically, and specialize at the non-CM Shimura
   point.

This description gives a reproducible reconstruction strategy, but not the
record-specific seed equations, finite-field solution, neighbor chain, or
final 17 sections. Recovering or independently rebuilding those objects is
the next geometry milestone.
