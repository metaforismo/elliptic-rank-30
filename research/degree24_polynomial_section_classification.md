# Complete classification of the degree-(2,4) coefficient scheme

**Status:** `proved` intermediate construction theorem.

This result classifies a restrictive polynomial identity.  It does not itself
produce an elliptic curve of rank 30.

## Classification theorem

Fix \(p\in\mathbf Q^*\).  Let

\[
L(v)=av^2+bv+c,\qquad a\ne0,
\]

and

\[
Q(v)=dv^4+ev^3+fv^2+gv+h.
\]

There exist rational coefficients and an \(S\in\mathbf Q\) satisfying

\[
Q(v)^2-v^2L(v)^3=v^3-Sv^2+p v+1 \tag{1}
\]

if and only if there is a \(t\in\mathbf Q^*\) such that

\[
\boxed{p=\Phi(t):=\frac{t^{12}}{186624}-\frac{t^3}{6}.} \tag{2}
\]

For such a parameter, every solution is obtained, after replacing \(Q\) by
\(-Q\) when necessary, from

\[
\begin{aligned}
a&=\frac{36}{t^4},
&b&=\frac{t^2}{6},
&c&=\frac{t^8}{5184}-\frac4t,\\
d&=\frac{216}{t^6},
&e&=\frac32,
&f&=\frac{t^6}{288}-\frac{36}{t^3},\\
g&=\frac{p}{2},
&h&=1,
&S&=c^3-g^2-2f.
\end{aligned} \tag{3}
\]

The sign that appears while taking a square root in the forward derivation is
absorbed by \(t\mapsto-t\), so (2)--(3) are a single complete family.

## Forward derivation

The constant and linear coefficients normalize to

\[
h=1,\qquad g=\frac p2.
\]

As in the exact verifier, the coefficients from degree eight down to degree
five give, for rational \(u,x\ne0\),

\[
\begin{aligned}
a&=u^2,&b&=ux,&d&=u^3,\\
e&=\frac32u^2x,
&c&=\frac{x^3+8p}{12x},
&f&=\frac{u(x^3+2p)}{2x}.
\end{aligned} \tag{4}
\]

Writing

\[
D=(x^3-4p)^2,
\]

the remaining degree-four and degree-three equations are exactly

\[
D=96ux^2, \tag{5}
\]

\[
144u^2x^2-uD-48x=0. \tag{6}
\]

The exact linear combination

\[
(6)-u(5)=48x(u^2x-1)
\]

forces

\[
u^2x=1. \tag{7}
\]

Equation (5) says that \(96u\) is a rational square.  Choose \(w\in\mathbf
Q^*\) with \(w^2=96u\), and absorb the sign of
\(x^3-4p=\pm wx\) into the parameter

\[
t=\pm\frac{24}{w}.
\]

Then

\[
u=\frac6{t^2},\qquad x=\frac{t^4}{36},
\]

and solving for \(p\) gives (2).  Substitution into (4) gives (3).

The assumption \(p\ne0\) excludes the boundary branch \(b=0\), which is why
\(x=b/u\) is legitimate.  The case needed for the historical boundary problem
is \(p=3\), so no relevant solution is lost.

## Converse

Substituting (2)--(3) into the left side of (1), with exact Laurent-polynomial
arithmetic over \(\mathbf Q[t,t^{-1},v]\), gives identically

\[
v^3-Sv^2+p v+1.
\]

The executable verifier performs this full symbolic substitution rather than
checking only numerical samples.

## Cube-coset representation

Put

\[
q=\frac{t^3}{36}.
\]

Then the classification map becomes

\[
\boxed{p=9q^4-6q,\qquad q\in\frac1{36}(\mathbf Q^*)^3.} \tag{8}
\]

This is the useful new representation: the coefficient problem is an image
problem for a quartic map restricted to one rational cube coset.  Equal values
of the quartic correspond to the cubic relation

\[
3(q_1+q_2)(q_1^2+q_2^2)=2, \tag{9}
\]

which is birational to the fixed Mordell curve

\[
y^2=x^3-1296. \tag{10}
\]

One convenient birational map is

\[
q_1=\frac{36+y}{6x},\qquad
q_2=\frac{36-y}{6x}.
\]

To produce two genuinely distinct decompositions for the same \(p\), a point
on (10) must additionally satisfy the two cube-coset conditions

\[
36q_1,36q_2\in(\mathbf Q^*)^3.
\]

Thus collision search is converted into an explicit elliptic-curve plus
3-descent problem rather than a blind coefficient search.

## Exact obstruction at the historical coefficient \(p=3\)

Equation (8) would require

\[
3q^4-2q-1=0.
\]

But

\[
3q^4-2q-1=(q-1)(3q^3+3q^2+3q+1).
\]

The cubic has no rational root by the rational-root theorem.  Hence the only
rational possibility is \(q=1\), which would require

\[
t^3=36.
\]

This is impossible in \(\mathbf Q\), since \(v_2(36)=2\not\equiv0\pmod3\).
Therefore the target coefficient \(p=3\) is not in the rational image of the
complete degree-(2,4) family.

## An exact sample

For \(t=3\),

\[
p=-\frac{423}{256},\qquad S=-\frac{665}{216},
\]

and

\[
L(v)=\frac49v^2+\frac32v-\frac{13}{192},
\]

\[
Q(v)=\frac8{27}v^4+\frac32v^3+\frac{115}{96}v^2
      -\frac{423}{512}v+1.
\]

Exact expansion gives

\[
Q(v)^2-v^2L(v)^3
=v^3+\frac{665}{216}v^2-\frac{423}{256}v+1.
\]

## Reproduction

```bash
python3 research/degree24_polynomial_section_classification.py \
  --check certificates/degree24_polynomial_section_classification.json
python3 -m unittest tests.test_degree24_polynomial_section_classification -v
```

The verifier checks the forward reduction with symbolic \(p\), the converse
identity with symbolic \(t\), five independent exact rational samples, the
factorization at \(p=3\), and the cube valuation obstruction.

## Research consequence

The degree-(2,4) boundary is not an unexplored high-dimensional coefficient
space.  It is a one-dimensional rational locus with an additional cube-class
restriction, and the historical value \(p=3\) lies outside that locus.  Future
inverse-construction work should therefore concentrate on:

1. collision points on (10) satisfying both cube conditions;
2. higher-degree polynomial identities;
3. rational functions with controlled poles;
4. fiber products of several distinct forcing maps;
5. descent structures that can certify independence immediately.
