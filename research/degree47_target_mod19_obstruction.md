# Pure mod-19 obstruction for the generic degree-(4,7) target

**Status:** `proved` exact obstruction.  No computer algebra system is needed
to verify the final theorem.

## Reduced system

In the generic Mason equality case, normalize

\[
Q(0)=1,\qquad [v]Q=\frac32,
\]

and write

\[
L=A(1+bv^2+cv^3+dv^4).
\]

The Wronskian identity determines \(Q\) recursively.  The coefficient of
\(v^3\) in the target identity forces

\[
c=\frac2{11}-3b.
\]

Introduce

\[
t=11b,\qquad z=121d,\qquad k=2-3t.
\]

The remaining Wronskian coefficients \(E_8,E_9,E_{10}\) satisfy the exact
identities

\[
11^4E_8=2H,\qquad
H=16z^2+198kz-33tk, \tag{1}
\]

\[
11^4E_9-9H=G_9, \tag{2}
\]

and

\[
11^6E_{10}-16t^2H=12kG_{10}, \tag{3}
\]

where

\[
G_9=(12t^2+5290t-3564)z
     -588t^3-4455t^2+6138t-1936, \tag{4}
\]

\[
G_{10}=(-264t^2+3993t-2904)z
       +44t^3-53845t^2+71874t-23958. \tag{5}
\]

If \(k=0\), equation (1) gives \(z=0\), hence \(d=0\), so \(L\) is not
quartic.  Thus a genuine solution has \(k
e0\), and (1)--(3) force

\[
G_9=G_{10}=0.
\]

## Determinant quintic

Eliminating \(z\) from (4)--(5) yields

\[
\begin{aligned}
P(t)={}&154704t^5-758384t^4+266432683t^3\\
      &-533872086t^2+357341556t-79764168.
\end{aligned} \tag{6}
\]

The executable certificate verifies (6) twice:

1. directly in a sparse Laurent-polynomial ring over \(\mathbf Q\);
2. by independent integer coefficient convolution of the two \(2\times2\)
   determinant products.

The polynomial \(P\) is primitive and

\[
154704\not\equiv0\pmod{19}.
\]

Its values at \(0,1,\ldots,18\) modulo \(19\) are

\[
17,10,8,11,1,16,15,9,2,8,11,1,10,15,12,14,11,6,13.
\]

None is zero.

If a primitive integer polynomial has a rational root \(m/n\) in lowest terms,
then \(n\) divides its leading coefficient.  Since \(19\) does not divide that
coefficient, \(n\) is invertible modulo \(19\), and \(m/n\) would reduce to a
root of \(P\) in \(\mathbf F_{19}\).  The exhaustive table proves that this is
impossible.

Therefore the generic degree-(4,7) target system has no rational solution.
Together with the separate repeated-root obstruction, there is no
\((4,7)\) polynomial solution at all.

## Global consequence within this representation

The complete degree-(2,4) classification excludes the target coefficient
\(p=3\).  Mason--Stothers proves that no levels beyond \((4,7)\) exist.  The
present theorem excludes the generic \((4,7)\) branch, while
`degree47_degenerate_obstruction.md` excludes the only repeated-root branches.
Hence:

\[
\boxed{\text{Every polynomial identity of this shape at }p=3\text{ is impossible.}}
\]

This closes the **polynomial** representation only.  It does not cover rational
functions with finite poles and does not imply any universal rank bound.

## Reproduction

```bash
python3 research/degree47_target_mod19_obstruction.py \
  --check certificates/degree47_target_mod19_obstruction.json
python3 -m unittest tests.test_degree47_target_mod19_obstruction -v
```
