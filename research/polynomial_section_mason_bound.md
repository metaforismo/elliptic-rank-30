# Mason--Stothers bound for the polynomial-section ansatz

**Status:** `proved` structural theorem.  This is not a rank-30 certificate.

## Theorem

Suppose

\[
Q(v)^2-v^2L(v)^3=R(v), \tag{1}
\]

where

\[
\deg L=2k,\qquad \deg Q=3k+1,\qquad \deg R\le3,
\]

and \(Q(0)R(0)\ne0\).  Then

\[
\boxed{k\le2.}
\]

Consequently, after the completely classified
\((\deg L,\deg Q)=(2,4)\) level, the **only** possible higher polynomial
level is

\[
\boxed{(\deg L,\deg Q)=(4,7).}
\]

No polynomial search of the same shape with \(\deg L\ge6\) can succeed.
Rational functions with poles are not covered by this theorem.

## Proof

Put

\[
A=Q^2,\qquad B=-v^2L^3,\qquad C=R,
\]

so that \(A+B=C\).  Let

\[
H=\gcd(A,B).
\]

Because \(H\mid A+B=R\), writing \(h=\deg H\) gives \(h\le3\).  Divide the
three terms by \(H\):

\[
A_1=\frac{Q^2}{H},\qquad
B_1=-\frac{v^2L^3}{H},\qquad
C_1=\frac RH.
\]

The polynomials \(A_1,B_1,C_1\) are pairwise coprime.  Their two large terms
have degree

\[
\deg A_1=\deg B_1=6k+2-h. \tag{2}
\]

The radical of \(A_1\) is supported on the roots of \(Q\), the radical of
\(B_1\) on the roots of \(vL\), and the radical of \(C_1\) on the roots of
\(R/H\).  Therefore

\[
\deg\operatorname{rad}(A_1B_1C_1)
 \le (3k+1)+(2k+1)+(3-h)
 =5k+5-h. \tag{3}
\]

Mason--Stothers applied to \(A_1+B_1=C_1\) yields

\[
6k+2-h
\le \deg\operatorname{rad}(A_1B_1C_1)-1
\le5k+4-h.
\]

Thus \(k\le2\).

## The maximal nondegenerate case as a Belyi problem

Take \(k=2\), so

\[
\deg L=4,\qquad \deg Q=7.
\]

Assume first that \(Q\), \(vL\), and \(R\) have disjoint squarefree supports.
Equality holds throughout the Mason estimate.  Define

\[
f(v)=\frac{Q(v)^2}{R(v)}.
\]

Then

\[
f(v)-1=\frac{v^2L(v)^3}{R(v)}.
\]

Hence the ramification partitions are:

- above \(0\): seven double zeros, \((2^7)\);
- above \(1\): four triple zeros and the double zero at \(v=0\),
  \((3^4,2)\);
- above \(\infty\): the pole of order \(11\) at infinity and three simple
  finite poles, \((11,1,1,1)\).

Thus the maximal polynomial problem is a finite genus-zero Belyi-map problem
with passport

\[
\boxed{(2^7),\quad(3^4\,2),\quad(11\,1^3).} \tag{4}
\]

The ramification contributions are

\[
7,\qquad9,\qquad10,
\]

whose sum \(26\) equals \(2\cdot14-2\), as required by Riemann--Hurwitz.

The same equality case gives a useful differential equation.  The Wronskian
of \(Q^2\) and \(-v^2L^3\) is divisible by \(vL^2Q\), of degree

\[
1+2\cdot4+7=16.
\]

Rewriting the Wronskian using the cubic remainder bounds its degree by

\[
14+3-1=16.
\]

Therefore the quotient is a nonzero constant, and

\[
\boxed{2LQ+3vL'Q-2vLQ'=\lambda,\qquad\lambda\in\mathbf Q^*.} \tag{5}
\]

Equation (5) is a far smaller exact system than the original coefficient
comparison.

## Degenerate target cubics

For the historical coefficient pattern

\[
R(v)=v^3-Sv^2+3v+1,
\]

the discriminant is

\[
\operatorname{Disc}(R)=4S^3+9S^2-54S-135
=(S+3)^2(4S-15). \tag{6}
\]

The only repeated-root parameters are

\[
S=-3,\qquad R=(v+1)^3,
\]

and

\[
S=\frac{15}{4},\qquad
R=(v-2)^2\left(v+\frac14\right).
\]

They must be checked separately because cancellation changes the generic
passport.  Every other \(S\) is governed by (4)--(5).

## Consequence for the search

The inverse-construction thread has become finite and structured:

1. enumerate the transitive permutation triples with passport (4);
2. determine their fields of moduli and fields of definition;
3. reconstruct the associated Belyi maps exactly;
4. test whether a rational normalization produces the coefficient \(3v\);
5. analyze the two degenerate values in (6) separately.

This is strictly stronger than increasing the degree of a blind polynomial
ansatz: Mason--Stothers proves that no later polynomial level exists.

## Reproduction

```bash
python3 research/polynomial_section_mason_bound.py \
  --check certificates/polynomial_section_mason_bound.json
python3 -m unittest tests.test_polynomial_section_mason_bound -v
```
