# Correction: the target-component F7 seed has no lift modulo 49

## Claim status

**Proved for one exact modular seed and one exact component/sign chart.**

The normalized semistable surface-section seed over \(\mathbf F_7\) does not
lift modulo \(7^2\) once the required split-\(I_4\) and split-\(I_3\) component
conditions are imposed as radical node-and-tangent equations. Therefore there
is no \(\mathbf Z_7\)-point in that chart reducing to this seed.

This does not exclude other modular seeds, other primes, the additive \(IV\)
realization, or a rank-30 curve.

## Why the previous nilpotence statement is withdrawn

The 19 coefficient equations consist of six equations for the surface and
thirteen coefficients of

\[
Y(t)^2=X(t)^3-3c_4(t)X(t)-2c_6(t).
\]

At the seed their Jacobian has rank 14. These equations record passage through
a repeated root via powers such as \(Y(0)^2\), rather than by the radical
conditions selecting the node and a tangent branch.

A previous note stated that the completed 19-equation local ideal contained

```text
b^3, w^4, t^12
```

and was therefore Artinian. That claim is withdrawn. The referenced generator
and certificate were never successfully generated, and an independent exact
truncation check does not support the asserted \(t^{12}\) membership.

The strategically relevant question is not whether the nonreduced coefficient
scheme has a long formal direction. It is whether that direction remains after
the required component data are imposed radical-theoretically.

## Radical component equations

For the constant denominator chart \(D=1\), the certified modular section has
tangent ratios \((+1,+1)\). The required local equations are

\[
\begin{aligned}
X(0)+p_0&=0, & X(1)+q_0&=0,\\
Y(0)&=0, & Y(1)&=0,\\
Y'(0)-X'(0)-p_1&=0,
&Y'(1)-X'(1)-q_1&=0.
\end{aligned}
\]

Together with the original 19 equations, these give 25 equations in 17
variables. Their Jacobian at the seed has full column rank 17 over
\(\mathbf F_7\).

## Exact mixed-characteristic obstruction

Suppose a lift modulo 49 existed. It would have the form

\[
x=x_0+7d \pmod {49},
\]

where \(x_0\) is the certified seed. Dividing the first-order expansion by 7
gives a linear system

\[
Jd=r \pmod 7.
\]

The executable certificate constructs an explicit vector \(\ell\) in the left
kernel of the \(25\times17\) Jacobian such that

\[
\ell J=0,
\qquad
\ell r=3\ne0 \pmod7.
\]

This is a contradiction. Hence

\[
\boxed{
\text{there is no lift modulo }49
\text{ of this seed in the required component chart.}
}
\]

A \(\mathbf Z_7\)-point reducing to the seed would reduce modulo 49, so it is
excluded as well.

## Reproduction

```bash
python3 research/rank17_component_lift_obstruction_f7.py \
  --compare certificates/rank17_component_lift_obstruction_f7.json

python3 -m unittest \
  tests.test_rank17_component_lift_obstruction_f7 -v
```

## Remaining boundary

This certificate closes only the exact semistable seed

```text
(p0,p1,p2,p3,q0,q1,q2,s,x0,x1,x2,x3,x4,y0,y1,y2,y3)
= (2,2,1,0,2,2,2,0,5,6,2,0,6,0,1,4,2)
```

with tangent ratios \((1,1)\) over \(\mathbf F_7\). It does not classify all
points of the semistable incidence locus. The additive \(IV\) branch and the
search for other characteristic-zero constructions remain open.
