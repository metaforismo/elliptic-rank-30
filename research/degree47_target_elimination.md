# Exact elimination of the generic degree-(4,7) target system

**Status:** executable exact computation.  The generated certificate decides
whether the generic maximal polynomial level contains a rational identity with
the historical coefficient \(3v\).

## Normalization

Assume the squarefree/coprime Mason equality case and write

\[
Q(v)^2-v^2L(v)^3=v^3-Sv^2+3v+1. \tag{1}
\]

Replace \(Q\) by \(-Q\) if necessary, so

\[
Q(0)=1,\qquad [v]Q=\frac32.
\]

The Wronskian equality is

\[
2LQ+3vL'Q-2vLQ'=\lambda. \tag{2}
\]

Its constant term is \(2L(0)\).  The coefficient of \(v\) forces the linear
coefficient of \(L\) to vanish.  Put

\[
L=A(1+bv^2+cv^3+dv^4),\qquad A\ne0. \tag{3}
\]

Equation (2) recursively determines all coefficients \(q_2,\ldots,q_7\) of
\(Q\) as rational functions of \(b,c,d\).  The coefficient of \(v^3\) in (1)
then gives

\[
\boxed{c=\frac2{11}-3b.} \tag{4}
\]

The coefficients of degrees \(8,9,10\) in (2) leave three exact equations in
only \(b,d\).  The program computes their lexicographic Gröbner basis over
\(\mathbf Q\), extracts and factors the elimination polynomial in \(b\), and
enumerates every rational common zero.

## Leading cube-class condition

For every rational pair \((b,d)\), the leading coefficient of \(Q\), denoted
\(q_7\), must satisfy

\[
q_7^2=(Ad)^3. \tag{5}
\]

Thus

\[
A^3=\frac{q_7^2}{d^3}. \tag{6}
\]

The verifier performs an exact numerator/denominator cube test.  Every
survivor is reconstructed as explicit \(L,Q,A,S\), and (1) is expanded over
\(\mathbf Q[v]\) coefficient by coefficient.  No numerical root finding enters
the certificate.

## Coverage

For a squarefree target cubic, any common zero of \(Q\) and \(vL\) would give a
square or cube factor in the cubic remainder, so the supports are disjoint.
Mason equality then forces squarefreeness and the Wronskian identity (2).
Therefore every generic rational degree-(4,7) target identity is represented in
the final `full_rational_solutions` array.

The repeated-root values \(S=-3\) and \(S=15/4\) are covered independently by
`certificates/degree47_degenerate_obstruction.json`.

## Reproduction

The symbolic engine is pinned in the integration workflow:

```bash
python3 -m pip install 'sympy==1.13.3'
python3 research/degree47_target_elimination.py \
  --output certificates/degree47_target_elimination.json
python3 research/degree47_target_elimination.py \
  --check certificates/degree47_target_elimination.json
```

The JSON records the exact Gröbner basis, elimination polynomial,
factorization, rational roots, candidate pairs, cube-class decisions, complete
reconstructions, SymPy version, and a canonical SHA-256 digest.

## Interpretation

- An empty `full_rational_solutions` list, together with the separate
  degenerate obstruction, closes every degree-(4,7) target identity.
- A nonempty list is an explicit exact polynomial construction and must be
  tested immediately for the resulting elliptic section's order and
  independence.
