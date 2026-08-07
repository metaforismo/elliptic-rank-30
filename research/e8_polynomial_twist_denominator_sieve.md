# A denominator barrier for visible twist sections on the split \(E_8\) surface

## Setting

Let \(E/\mathbb Q(t)\) be the explicitly certified Kumar--Shioda rational
elliptic surface of Mordell--Weil rank eight, written

\[
E:\quad y^2=x^3+t^2x^2+a_4(t)x+a_6(t).
\]

For a branch pair \(\{a,b\}\subset\{-4,-3,\ldots,4,\infty\}\), let
\(d_{a,b}(t)\) be the square-free quadratic character ramified at that pair.
We ask for the lowest natural-degree polynomial sections of the twist

\[
d_{a,b}(t)y(t)^2=x(t)^3+t^2x(t)^2+a_4(t)x(t)+a_6(t),
\]

with

\[
\deg x\le2,
\qquad
\deg y\le2.
\]

This ansatz is the first place where a new low-height direction could become
visible without introducing poles.

## Exact exhaustive computation

For each prime

\[
23,29,31,37,41,47,53,61,71,73,79,
\]

the certificate enumerates all \(p^3\) possible triples
\((A,B,C)\in\mathbb F_p^3\) for

\[
x(t)=At^2+Bt+C.
\]

For each of the 45 branch pairs it then:

1. computes the exact right-hand polynomial;
2. checks divisibility by \(d_{a,b}(t)\);
3. divides exactly when possible;
4. tests whether the quotient is a polynomial square of degree at most four;
5. records every surviving \(x\)-class.

No random sampling, coefficient cutoff, floating point, or rank heuristic is
used.

The first four primes exhibit only sporadic residue-field sections:

\[
\begin{array}{c|c|c}
p&\text{nonempty branch pairs}&\text{total }x\text{-classes}\\
\hline
23&1&1\\
29&1&4\\
31&4&6\\
37&8&9.
\end{array}
\]

At each of

\[
41,47,53,61,71,73,79
\]

all 45 branch pairs are empty.

## Theorem

For every branch pair in

\[
\{-4,-3,\ldots,4,\infty\},
\]

any section over \(\mathbb Q(t)\) satisfying the polynomial degree bounds
\(\deg x\le2\), \(\deg y\le2\) must have a common coefficient denominator
divisible by

\[
41\cdot47\cdot53\cdot61\cdot71\cdot73\cdot79
=
\boxed{2550913424887}.
\]

### Proof

Suppose such a rational section had coefficients integral at one of the seven
listed primes. The Kumar--Shioda model has good coefficient reduction there.
Reducing the section coefficientwise would give a solution of the same
polynomial identity over \(\mathbb F_p(t)\) in the exhaustively enumerated
ansatz. The exact computation proves that no such solution exists. Therefore
the section must have a coefficient denominator divisible by that prime.
Applying the argument independently at all seven primes gives the stated
product divisor.

## Interpretation

The finite-field sections visible at small primes are not stable arithmetic
signals; they disappear completely at seven other good primes. The failure is
therefore not merely that a coefficient box was too small. On this explicit
split \(E_8\) surface, every section in the lowest polynomial-height shell is
forced to have an enormous denominator.

This explains why a surface may have excellent lattice capacity and still be
poor for rational-point visibility. The next search must move to one of:

- rational-function sections with controlled poles;
- a different point of the split \(E_8\) moduli family with smaller arithmetic
  coefficients;
- branch supports selected by a denominator-aware objective;
- covering curves or descent coordinates in which the same sections have lower
  height.

## Scope

This theorem is deliberately restricted. It does **not** prove that any twist
has rank zero, and it does not exclude sections of higher polynomial degree or
sections with poles. It is an exact visibility obstruction for one explicit
surface, one branch-support box, and the natural minimal polynomial ansatz.

The executable computation is
`search/e8_polynomial_twist_denominator_sieve.py`; the compact certificate is
`certificates/e8_polynomial_twist_denominator_sieve.json`.
