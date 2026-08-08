# Current record mechanism: verified facts, inherited search machinery, and missing reconstruction

## Scope

This document records only facts supported by primary sources or by exact
certificates in this repository. It deliberately separates:

- the **unconditional rank-29 lower bound**;
- the **conditional exact-rank statement**;
- the publicly described search mechanism;
- the record-specific geometric data that have not yet been recovered.

Primary sources:

1. Noam Elkies, “\( \mathbf Z^{29} \) in \(E(\mathbf Q)\),” Number Theory
   List, 29 August 2024:
   <https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY%3Bb9d018b1.2409&S=>
2. Noam D. Elkies and Zev Klagsbrun, *New Rank Records for Elliptic Curves
   Having Rational Torsion*, arXiv:2003.00077:
   <https://arxiv.org/pdf/2003.00077>
3. Zev Klagsbrun, Travis Sherman, and James Weigandt, *The Elkies
   Rank Bound Algorithm*, Math. Comp. 88 (2019), 837–846:
   <https://arxiv.org/abs/1606.07178>

## 1. Record curve and proof status

The published integral model is

\[
E_{29}:\quad
y^2+xy
=
x^3
-27006183241630922218434652145297453784768054621836357954737385\,x
\]
\[
\phantom{E_{29}:\quad y^2+xy={}}
+55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497.
\]

The 2024 announcement states:

- \(E_{29}(\mathbf Q)\) contains a subgroup of rank 29;
- the curve has rank exactly 29 assuming GRH for zeta functions of number
  fields;
- an analytic-rank upper bound of 29 is also conditional on an \(L\)-function
  GRH, with BSD needed to identify analytic and algebraic rank in that route.

Only the lower bound is needed for the record and for this project’s
baseline.

### Repository verification

The files

- `curve.json`;
- `points.json`;
- `baseline/verify_rank29_mod2.py`;
- `baseline/rank29_mod2_certificate.json`;

give a new exact standard-library verification of

\[
\operatorname{rank}E_{29}(\mathbf Q)\ge29.
\]

The proof checks all 29 points by rational substitution, proves
\(E_{29}(\mathbf Q)_{\rm tors}=0\) from reductions at 67 and 71, and verifies
that the 29 local images have rank 29 in a product of
\(E(\mathbf F_p)/2E(\mathbf F_p)\). It uses no floating point arithmetic,
GRH, BSD, SageMath, Magma, or PARI/GP.

This does not reproduce the conditional upper bound and does not claim that
the 29-point subgroup is globally saturated.

## 2. Record-specific geometry known publicly

The public announcement gives the following geometric accounting:

- the search used a **rank-17 elliptic fibration** of the same K3 surface that
  produced the 2006 rank-28 record;
- for each high-scoring specialization, points were sought outside the
  generic \(\mathbf Z^{17}\);
- the record specialization produced **12 additional independent points**.

Thus the observed lower-bound decomposition is

\[
29=17+12.
\]

The target rank 30 requires one of two broad improvements:

\[
17+13,
\]

or a different fibration/construction with a larger generic floor and a
smaller exceptional jump.

The announcement says that the detailed K3 computation and search write-up
were intended for later publication. The explicit rank-17 fibration,
specialization parameter, 17 generic sections, and specialization map are
not contained in the announcement itself. They must not be invented or
inferred from the final curve alone.

## 3. Search machinery described by Elkies--Klagsbrun

The 2020 paper describes the search architecture later cited in the rank-29
announcement.

### 3.1 Generic rank as a floor

Start with an elliptic fibration

\[
\mathcal E/\mathbf Q(t)
\]

of Mordell--Weil rank \(r\). For all but finitely many specializations,
Silverman’s specialization theorem preserves at least those \(r\) directions.
The specialization search therefore hunts for an exceptional jump above a
known generic subgroup instead of starting from rank zero.

### 3.2 Mestre--Nagao score

For a specialization \(t\) and prime bound \(B\), the paper uses

\[
S(t,B)
=
\sum_{\substack{p<B\\E_t\ {\rm good\ at}\ p}}
\log\!\left(\frac{\#E_t(\mathbf F_p)}{p}\right).
\]

This is a ranking heuristic only. A large score triggers more expensive
arithmetic; it is never a rank certificate.

### 3.3 Precompute local traces

For each prime \(p\), the trace \(a_p(E_t)\) depends only on \(t\bmod p\).
The scanner can therefore precompute all local contributions for every
residue class before evaluating a large rational search region.

### 3.4 Staged cutoffs

Choose increasing prime bounds

\[
B_0\le B_1\le\cdots\le B_m=B
\]

and thresholds

\[
C_0\le C_1\le\cdots\le C_m.
\]

A parameter survives to stage \(i\) only if it passed every earlier score
threshold. This prevents expensive large-\(B\) evaluations on the full
search region.

### 3.5 Sieve over rational parameters

For fixed denominator \(b\) and a long interval of numerators \(a\), the
score contributions are added to an array for all \(t=a/b\) at once.
The paper reports:

- fixed-point approximations to
  \(\log(\#E_t(\mathbf F_p)/p)\);
- denominator \(D=1024\);
- 16-bit score counters;
- SIMD/vector additions;
- the initial sieve as the dominant bulk stage.

These are implementation choices to reproduce and benchmark, not immutable
constants.

### 3.6 Descent before deep point search

After scoring, the 2020 pipeline applies descent/Selmer computations to
obtain stronger rank information or rejection filters where practical. The
exact descent depends on available torsion/isogeny structure.

For survivors, rational points are searched on covering curves, especially
2-coverings, rather than only by direct bounded \(x\)-coordinate searches.
When known points give new coverings, those coverings are added to the search.

### 3.7 Skewed parameter regions

If the coefficient polynomials are unbalanced, square boxes in numerator and
denominator need not minimize the size of specialized coefficients. The paper
uses skewed regions to favor smaller resulting curves at comparable parameter
height.

## 4. What likely mattered structurally

The following are research hypotheses, not facts about the record until they
are tested against the recovered fibration.

### H1. The score was identifying an exceptional local signature, not merely a
large average.

The rank-29 parameter should be compared with controls prime by prime.
Ablations must determine which primes contribute genuine predictive
information for actual extra points.

### H2. Bad-reduction data may carry signal omitted by the classical score.

The 2020 paper explicitly notes that split multiplicative reduction and
Tamagawa behavior are not naturally represented in \(S(t,B)\), and reports
that ad hoc bonuses sometimes helped. The record curve has a highly
structured discriminant and many multiplicative primes; whether this
structure predicted the twelve extra points is testable once the family is
recovered.

### H3. The fibration may contain hidden character packets.

The verified theorem in `research/genus_zero_three_channel_packet.md` shows
that two quadratic multisections can hide a third independent product-
character section on a genus-zero biquadratic base. Every pair of quadratic
multisections on the record fibration must therefore be audited across all
three nontrivial character twists.

### H4. Coefficient size and point height must be separated from score.

A score can remain large while the expected height of new points becomes too
large for the covering search. Search regions should be evaluated by cost per
additional independent point, not score alone.

## 5. Reconstruction acceptance criteria

The current-record geometry is considered reproduced only when `main`
contains exact data and scripts satisfying all of the following.

### Geometry

- an explicit Weierstrass model over \(\mathbf Q(t)\);
- the K3 surface/pencil from which it is derived;
- singular fibers and Shioda--Tate accounting;
- 17 explicit sections over \(\mathbf Q(t)\);
- an exact height-pairing matrix of rank 17;
- the specialization parameter producing \(E_{29}\);
- the exact isomorphism from the specialized fiber to `curve.json`.

### Point accounting

- specialization of all 17 generic sections;
- identification of 12 additional independent directions;
- an exact transition matrix from those 29 directions to the published
  29-point basis;
- checks that no point is introduced over a number field instead of
  \(\mathbf Q\).

### Search reproduction

- the local-score table for every prime used;
- exact \(B_i,C_i\), fixed-point scale, and search region;
- a deterministic scan that ranks the known record parameter;
- point-search/covering commands and seeds;
- ablations connecting score features to actual extra points.

### Independent verification

- the existing exact mod-2 certificate;
- a SageMath verification script;
- a Magma verification script;
- software versions and logs;
- a clear separation between unconditional lower bounds and conditional
  upper bounds.

## 6. Immediate executable experiments after recovery

1. **Prime ablation:** remove one prime or one prime band at a time and measure
   the change in record-parameter rank among controls.
2. **Bad-prime features:** add only interpretable reduction-type/Tamagawa
   features and test out-of-sample point yield.
3. **Character audit:** enumerate squareclass spans of quadratic
   multisections; test every product twist.
4. **Exceptional-point height model:** regress actual covering-search cost and
   point height against local features, not analytic-rank labels.
5. **Neighbor fibrations:** reconstruct the Néron--Severi lattice and
   enumerate alternative \(U\)-embeddings/neighbor steps, ranking them by
   generic rank, rationality of sections, and coefficient growth.

## 7. Current blocker

The exact lower-bound verification is complete. The main blocker is no longer
the final curve or its 29 points; it is the missing bridge back to the
rank-17 K3 fibration:

\[
\boxed{
\text{explicit fibration}
+\text{17 sections}
+\text{record parameter}
+\text{search provenance}.
}
\]

Until that bridge is recovered from primary material or independently
reconstructed, a production scan claiming to extend the original search
would not be reproducible.
