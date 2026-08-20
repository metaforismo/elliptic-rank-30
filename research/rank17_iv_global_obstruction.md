# Global obstruction for the normalized additive-\(IV\) branch

```text
SOLVED-A: false
SOLVED-B: false
unconditional global lower bound: rank E(Q) >= 29
```

## Result

There is no nondegenerate elliptic K3 surface over \(\mathbf Q\) in the
normalized

\[
I_{12}+I_4+IV
\]

locus for which the \(I_4\) and \(IV\) fibres are simultaneously split in
the classes required by the rank-17 \(X(6,79)\) seed lattice.

This eliminates the additive-\(IV\) realization of the \(A_2\) fibre from
the current rank-30 strategy. It does **not** eliminate the semistable
\(I_3\) realization.

## Algebraic component classification

The logarithmic-derivative reduction has finite projection to

\[
(k,r_3)
\]

because the first residual equation \(e_8\) has degree \(4\) in \(r_2\)
with constant leading coefficient

\[
23808556800.
\]

The exact common curve factor of

\[
\operatorname{Res}_{r_2}(e_8,e_7),\qquad
\operatorname{Res}_{r_2}(e_8,e_6),\qquad
\operatorname{Res}_{r_2}(e_8,e_5)
\]

is

\[
\boxed{F\,S^2},
\]

where \(F,S\in\mathbf Q[k,r_3]\) are irreducible quintics.

A degree-one subresultant identifies the generic lift over \(F\) with the
known rational component \(R_1=R_2=0\). The same subresultant identifies
the generic lift over \(S\) with the determinant-zero equation

\[
D=
-2k^2-37kr_3+52k+108r_2
-116r_3^2+148r_3-320=0.
\]

Since the projection is finite, every positive-dimensional component has
positive-dimensional image. Hence every such component lies over \(F\)
or \(S\).

## The determinant-zero component has no rational point

On \(D=0\),

\[
r_2=
\frac{
2k^2+37kr_3-52k+116r_3^2-148r_3+320
}{108}.
\]

The next two differential coefficients determine \(r_1,r_0\), while the
following coefficient is exactly

\[
\frac7{135}S(k,r_3).
\]

Modulo \(S\), the next equation is

\[
11907\,L(k,r_3)m_0+C(k,r_3)=0.
\]

The special chart \(L=0\) is empty, because

\[
\gcd\left(
\operatorname{Res}_k(S,L),
\operatorname{Res}_k(S,C)
\right)=1.
\]

On \(L\ne0\), solving for \(m_0\) and imposing the next coefficient gives
a polynomial \(G(k,r_3)\). The exact resultant

\[
\operatorname{Res}_k(S,G)
\]

has two irreducible factors over \(\mathbf Q\), of degrees

\[
9,\qquad56.
\]

It has no linear factor. Thus no rational value of \(r_3\), and hence no
rational affine point \((k,r_3)\), can occur on the determinant-zero
chart.

## The rational component cannot be simultaneously split

The other component is rational, with parameter \(w\). Its two split
targets have square classes

\[
\frac{w+3}{(w-3)(w^2+3)}
\]

and

\[
-\frac{2w(w+3)}{w^2+3}.
\]

A nonzero rational square is positive over \(\mathbf R\). The first class
is positive only on

\[
(-\infty,-3)\cup(3,\infty),
\]

whereas the second is positive only on

\[
(-3,0).
\]

The intervals are disjoint. All exceptional parameters are boundary or
were checked separately.

Therefore \(F\) contains no nondegenerate rational point satisfying both
split conditions.

## Why no isolated rational surface is missed

The nondegenerate surface defines a degree-\(20\) cover with ordered
passport

\[
(3^6,2),\quad
(2^{10}),\quad
(12,4,1^4),\quad
(2,1^{18}).
\]

After fixing the first three branch values at \(0,1,\infty\), its ordered
Hurwitz space is an unramified cover of

\[
\mathbf P^1\setminus\{0,1,\infty\}.
\]

The uniquely ramified points of indices \(2,4,12\) rigidify the source
coordinate. Hence the normalized coefficient locus is locally a
one-dimensional Hurwitz chart and has no isolated nondegenerate
geometric points.

The full argument is in:

```text
research/rank17_iv_hurwitz_purity_lemma.md
```

Thus every rational open point must lie on one of the
positive-dimensional components already classified.

## Final consequence

Combining the three exact statements:

1. every nondegenerate rational point lies on a curve component;
2. every curve component lies over \(F\) or \(S\);
3. \(F\) fails the split conditions and \(S\) has no rational point;

gives

\[
\boxed{
U_{I_{12}+I_4+IV}^{\mathrm{split}}(\mathbf Q)=\varnothing.
}
\]

The next construction-first target is the semistable

\[
I_{12}+I_4+I_3
\]

branch, where the surface Hurwitz locus is two-dimensional and the
height-\(79/12\) section should cut out a one-dimensional incidence locus.

## Evidence files

- `research/data/RANK17_IV_COMPONENT_CLASSIFICATION_FROZEN.json`
- `research/verify_rank17_iv_component_classification.py`
- `research/rank17_iv_rational_component_split_obstruction.py`
- `research/rank17_iv_three_variable_elimination_sage.py`

The frozen verifier has a fast mode and an optional exact heavy-resultant
mode.
