# Exact Lagrange-resolvent hash for cyclic cubic covers

**Author:** Francesco Giannicola  
**Truth status:** exact algebraic reduction; no rank-30 curve is claimed.

## Depressed cubic

Let `K` be a field of characteristic different from two and three. Start from
an irreducible monic cubic

\[
f(x)=x^3+A x^2+B x+C.
\]

After `x=X-A/3`, write

\[
f(X-A/3)=X^3+pX+q,
\]

where

\[
p=B-\frac{A^2}{3},
\qquad
q=\frac{2A^3}{27}-\frac{AB}{3}+C.
\]

Its discriminant is

\[
\Delta=-4p^3-27q^2.
\]

Assume `Delta=s^2` in `K`. Then the cubic has Galois group contained in `A3`;
when it is irreducible, its splitting field over `K` is cyclic of degree
three.

## Lagrange resolvents

Let the roots be `r1,r2,r3`, and let `zeta` be a primitive cube root of unity.
Put

\[
R=r_1+\zeta r_2+\zeta^2 r_3,
\qquad
S=r_1+\zeta^2 r_2+\zeta r_3.
\]

Then

\[
RS=-3p,
\qquad
R^3+S^3=-27q,
\]

and

\[
\boxed{
R^3=\frac{-27q+3s\sqrt{-3}}{2},
\qquad
S^3=\frac{-27q-3s\sqrt{-3}}{2}.
}
\]

The two right-hand sides are conjugate and their product is

\[
(-3p)^3=-27p^3.
\]

Therefore, over `K(zeta3)`, the cyclic cubic extension is the Kummer extension

\[
\boxed{
K(\zeta_3)\left(\sqrt[3]{d}\right),
\qquad
d=\frac{-27q+3s\sqrt{-3}}{2}.
}
\]

Changing the sign of `s` replaces `d` by its conjugate; modulo cubes this
corresponds to the inverse character. Thus the unordered pair of Kummer
classes

\[
\boxed{\{[d],[d]^{-1}\}}
\]

is a canonical character hash for the cyclic cubic extension.

## Application to trisection packets

For every irreducible square-discriminant cubic produced by a trace class:

1. depress the cubic exactly;
2. choose an exact square root `s` of the discriminant;
3. compute the Lagrange class `d` over `K(zeta3)`;
4. remove cube factors from its divisor;
5. store the unordered pair `[d],[d]^{-1}`.

Two trisections over the same base coordinate can define the same nontrivial
character channel only if their normalized Lagrange hashes agree. This is a
far cheaper first grouping test than constructing an isomorphism between the
full genus-one covers.

When the normalized divisor of `d` has three branch points and the extension
genus is one, the cover is a cubic twist of a `j=0` elliptic curve. Rational
descent and positive Mordell--Weil rank remain separate arithmetic tests.

The symbolic identities are verified in
`research/cyclic_cubic_lagrange_resolvent_certificate.py`.
