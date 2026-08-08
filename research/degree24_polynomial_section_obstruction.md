# Degree-(2,4) polynomial-section obstruction

**Status:** `proved` intermediate negative result.  This document does not claim
rank 30 and does not constrain arbitrary elliptic curves or arbitrary sections.

## Theorem

Let

\[
L(v)=av^2+bv+c,\qquad a\ne0,
\]

and

\[
Q(v)=dv^4+ev^3+fv^2+gv+h
\]

have coefficients in \(\mathbf Q\).  There is no \(S\in\mathbf Q\) for which

\[
Q(v)^2-v^2L(v)^3=v^3-Sv^2+3v+1. \tag{1}
\]

Thus the restrictive polynomial-section ansatz encoded by (1), with
\((\deg L,\deg Q)=(2,4)\), has no rational solution.

## Exact proof

The constant and linear coefficients of (1) give

\[
h^2=1,\qquad 2gh=3.
\]

Replacing \(Q\) by \(-Q\) if necessary, take

\[
h=1,\qquad g=\frac32.
\]

The coefficients of \(v^8\) through \(v^5\) are

\[
\begin{aligned}
d^2&=a^3,\\
2de&=3a^2b,\\
e^2+2df&=3a^2c+3ab^2,\\
2dg+2ef&=6abc+b^3.
\end{aligned}
\]

Because \(a\ne0\), also \(d\ne0\).  Put

\[
u=\frac da.
\]

Then \(a=u^2\) and \(d=u^3\).  The next two equations yield

\[
e=\frac32ub,
\qquad
f=\frac32uc+\frac38\frac{b^2}{u}.
\]

The coefficient of \(v^5\) becomes

\[
b^3-12u^2bc+24u^3=0.
\]

It cannot have \(b=0\), so write \(b=ux\) with \(x\ne0\).  Hence

\[
c=\frac{x^3+24}{12x},\qquad
 e=\frac32u^2x,\qquad
 f=\frac{u(x^3+6)}{2x}. \tag{2}
\]

Let

\[
D=(x^3-12)^2.
\]

After substituting (2), the coefficients of \(v^4\) and \(v^3\) reduce exactly
to

\[
D=96ux^2, \tag{3}
\]

and

\[
144u^2x^2-uD-48x=0. \tag{4}
\]

Subtracting \(u\) times (3) from (4) gives

\[
48x(u^2x-1)=0.
\]

Since \(x\ne0\),

\[
x=\frac1{u^2}.
\]

Substitution into (3), followed by setting \(z=u^3\), gives

\[
(1-12z^2)^2=96z^3,
\]

or equivalently

\[
144z^4-96z^3-24z^2+1=0. \tag{5}
\]

The exact factorization is

\[
144z^4-96z^3-24z^2+1
=(6z-1)(24z^3-12z^2-6z-1). \tag{6}
\]

The cubic factor has no root modulo \(7\): its values at
\(0,1,\ldots,6\) are

\[
6,5,5,3,3,2,4.
\]

Therefore it has no rational root.  The only remaining rational root of (5)
is \(z=1/6\), but this is not a cube in \(\mathbf Q\), since

\[
v_2(1/6)=-1\not\equiv0\pmod3.
\]

This contradicts \(z=u^3\), proving the theorem.

The coefficient of \(v^2\) merely determines

\[
S=c^3-\frac94-2f,
\]

so it introduces no omitted branch.

## Reproduction

From the repository root:

```bash
python3 research/degree24_polynomial_section_obstruction.py \
  --check certificates/degree24_polynomial_section_obstruction.json
python3 -m unittest tests.test_degree24_polynomial_section_obstruction -v
```

The verifier independently reconstructs all nine coefficient equations using a
tiny sparse Laurent-polynomial ring over \(\mathbf Q\), checks the substitutions
and elimination identities, verifies the integer factorization, enumerates the
cubic modulo \(7\), and checks the valuation obstruction.

## Consequence for the research program

This closes one tempting boundary representation: a rational section cannot be
obtained from identity (1) with quadratic \(L\) and quartic \(Q\).  The result
should redirect the inverse-construction thread toward at least one of:

1. a higher-degree \((L,Q)\) ansatz;
2. rational functions with controlled poles rather than polynomials;
3. a different norm form;
4. a multivariable auxiliary surface whose specialization produces sections.

It does **not** refute a broader modular or isogeny construction, and it does
not imply an upper bound on Mordell--Weil rank.
