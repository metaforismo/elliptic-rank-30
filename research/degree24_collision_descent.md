# Collision descent to simultaneous Gaussian and Eisenstein ninth powers

**Status:** `proved` necessary-condition theorem.  Nonexistence of the final
simultaneous ninth-power system remains open in this workstream.

## From a quartic collision to an integral equation

The complete degree-(2,4) classification uses

\[
q=\frac{t^3}{36},\qquad \Phi(t)=9q^4-6q.
\]

For distinct parameters \(t_1,t_2\), the equality
\(\Phi(t_1)=\Phi(t_2)\) is equivalent to

\[
3(q_1+q_2)(q_1^2+q_2^2)=2.
\]

Substituting \(q_i=t_i^3/36\) gives

\[
(t_1^3+t_2^3)(t_1^6+t_2^6)=31104=2^7 3^5. \tag{1}
\]

After clearing a common denominator and dividing common factors, any rational
collision gives primitive integers \(a,b,c\), with \(\gcd(a,b)=1\), satisfying

\[
(a^3+b^3)(a^6+b^6)=2^7 3^5 c^9. \tag{2}
\]

## Cyclotomic factors

Set

\[
\begin{aligned}
A&=a+b,\\
B&=a^2-ab+b^2,\\
C&=a^2+b^2,\\
D&=a^4-a^2b^2+b^4.
\end{aligned}
\]

Then

\[
a^3+b^3=AB,\qquad a^6+b^6=CD. \tag{3}
\]

The following exact identities control all common divisors:

\[
\begin{aligned}
B-(a+b)(a-2b)&=3b^2,\\
C-(a+b)(a-b)&=2b^2,\\
D-(a+b)a^2(a-b)&=b^4,\\
C-B&=ab,\\
D-(a^2+ab-b^2)B&=2b^3(b-a),\\
D-(a^2-2b^2)C&=3b^4.
\end{aligned} \tag{4}
\]

Using \(\gcd(a,b)=1\), these imply that the four factors are pairwise
coprime outside the indicated primes:

- \(\gcd(A,B)\mid3\);
- \(\gcd(A,C)\mid2\);
- \(\gcd(A,D)=\gcd(B,C)=\gcd(B,D)=\gcd(C,D)=1\)
  under the local conditions below.

## Forced local conditions

If one of \(a,b\) were even and the other odd, the left side of (2) would be
odd, contradicting its 2-adic valuation.  Hence \(a,b\) are both odd.  Complete
enumeration modulo \(8\) gives

\[
v_2(C)=1,\qquad v_2(B)=v_2(D)=0.
\]

Therefore

\[
v_2(A)=6+9v_2(c). \tag{5}
\]

Similarly, primitivity and (2) force \(3\nmid ab\) and
\(a+b\equiv0\pmod3\).  Enumeration modulo \(9\) gives

\[
v_3(B)=1,\qquad v_3(C)=v_3(D)=0,
\]

and consequently

\[
v_3(A)=4+9v_3(c). \tag{6}
\]

Every prime \(\ell\ge5\) occurs in exactly one of \(A,B,C,D\), and its
exponent there must be divisible by nine.

## Ninth-power system

Equations (5)--(6), pairwise coprimality, positivity of \(B,C,D\), and the
constant \(2^7 3^5\) in (2) force the simultaneous system

\[
\boxed{a+b=2^6 3^4 w^9,} \tag{7}
\]

\[
\boxed{a^2-ab+b^2=3x^9,} \tag{8}
\]

\[
\boxed{a^2+b^2=2y^9,} \tag{9}
\]

\[
\boxed{a^4-a^2b^2+b^4=z^9.} \tag{10}
\]

Conversely, multiplying (7)--(10) recovers the prime allocation in (2).
Thus every nontrivial collision must pass this much narrower Diophantine gate.

## Foreign-method interpretation

Equation (8) is an Eisenstein norm:

\[
a^2-ab+b^2=N_{\mathbf Q(\sqrt{-3})/\mathbf Q}(a+b\omega).
\]

Equation (9) is a Gaussian norm:

\[
a^2+b^2=N_{\mathbf Q(i)/\mathbf Q}(a+ib).
\]

Since both rings of integers are unique-factorization domains, primitive
solutions would force, after extracting the ramified factors above, simultaneous
ninth powers in \(\mathbf Z[\omega]\) and \(\mathbf Z[i]\).  The collision
problem has therefore become an intersection problem between two explicitly
parameterized generalized-Fermat loci rather than a blind rational search.

This opens several rigorous next routes:

1. derive the Gaussian ninth-power binary forms from
   \((a+ib)/(1+i)=\varepsilon(r+is)^9\);
2. derive the Eisenstein ninth-power forms from division by \(1-\omega\);
3. eliminate \(a,b\) between the two parameterizations;
4. use Thue--Mahler solvers or a modular/Frey-curve argument;
5. search for a single local obstruction on the resulting covering curves.

## Reproduction

```bash
python3 research/degree24_collision_descent.py \
  --check certificates/degree24_collision_descent.json
python3 -m unittest tests.test_degree24_collision_descent -v
```

The verifier checks the identities (3)--(4) symbolically and exhausts all
relevant residue classes modulo \(8\) and \(9\).

## Scope

This proves a necessary condition for two distinct degree-(2,4)
decompositions.  It does not yet prove that the ninth-power system has no
solutions, nor does it imply a global rank bound.
