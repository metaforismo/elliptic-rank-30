# Current record mechanism: verified facts, reconstructed public machinery, and missing bridge

## Scope

This document separates four layers that must not be conflated:

1. the **unconditional rank-29 lower bound** for the published curve;
2. the **conditional exact-rank statement**;
3. the rank-17 K3/ Shimura geometry and search machinery described in primary
   sources;
4. the record-specific equation, sections, parameter, and transformation data
   that are still missing.

Primary sources:

1. Noam Elkies, “\(\mathbf Z^{29}\) in \(E(\mathbf Q)\),” Number Theory
   List, 29 August 2024:
   <https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY%3Bb9d018b1.2409&S=>
2. Noam D. Elkies, *Three lectures on elliptic surfaces and curves of high
   rank*, arXiv:0709.2908:
   <https://arxiv.org/abs/0709.2908>
3. Noam D. Elkies, *Shimura curve computations via K3 surfaces of
   Néron-Severi rank at least 19*, arXiv:0802.1301:
   <https://arxiv.org/abs/0802.1301>
4. Noam D. Elkies and Zev Klagsbrun, *New Rank Records for Elliptic Curves
   Having Rational Torsion*, arXiv:2003.00077:
   <https://arxiv.org/pdf/2003.00077>
5. Zev Klagsbrun, Travis Sherman, and James Weigandt, *The Elkies Rank Bound
   Algorithm*, Math. Comp. 88 (2019), 837-846:
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

Only the lower bound is needed for the record and for this project's
baseline.

### Repository verification

The files

- `curve.json`;
- `points.json`;
- `baseline/verify_rank29_mod2.py`;
- `baseline/rank29_mod2_certificate.json`;

give an exact standard-library proof of

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

## 2. Geometric source of the rank-17 floor

### 2.1 Why rank 17 is the K3 target

For an elliptic K3 surface, the Néron-Severi rank is at most 20. The zero
section and fiber contribute a hyperbolic plane \(U\), while reducible fibers
consume root-lattice rank. Thus an elliptic K3 surface over \(\mathbf Q(t)\)
with no reducible fibers can have Mordell-Weil rank at most 18.

Elkies explains that rank 18 cannot occur over \(\mathbf Q\) in the desired
form, while rank 17 is attainable through a Néron-Severi lattice of rank 19
whose essential lattice is rootless. The corresponding moduli space is a
curve rather than an isolated CM point.

### 2.2 The exact Shimura anchor

The successful moduli curve has level

\[
N=6\cdot79=474
\]

and is the genus-two Shimura quotient

\[
X(6,79)/\langle w_{474}\rangle.
\]

Its published equation is

\[
C:\quad u^2=16t^6-19t^4+88t^2-48.
\]

The source lists two rational points at infinity, four affine points with
\(|t|=2\), \(|u|=32\), and four with

\[
|t|=\frac{14}{13},
\qquad
|u|=\frac{2^6\cdot251}{13^3}
=
\frac{16064}{2197}.
\]

It states that the last orbit is non-CM and yields an elliptic K3 surface of
Mordell-Weil rank 17 over \(\mathbf Q(s)\).

The repository now independently certifies the algebraic anchor in

- `research/rank17_shimura_anchor.py`;
- `certificates/rank17_shimura_anchor.json`;
- `research/rank17_shimura_anchor.md`.

The certificate proves that the sextic is squarefree and genus two, verifies
all eight affine rational points and the two rational points at infinity, and
proves that

\[
(t,u)\longmapsto(-t,-u)
\]

has the nonsingular genus-one quotient

\[
y^2=16x^4-19x^3+88x^2-48x,
\qquad x=t^2,\quad y=tu.
\]

It deliberately does **not** prove the non-CM status or reconstruct the K3
moduli map.

### 2.3 Lattice and p-adic construction route

The 2007 lectures describe how the original K3 equation and its generators
were computed.

1. **Essential-lattice design.** A suitable positive-definite essential
   lattice is obtained as a slice of a Niemeier lattice. Roots and glue control
   reducible fibers, Mordell-Weil torsion, and the target height pairing.
2. **Finite-field seed.** The polynomial identities for the elliptic surface
   and sections are solved modulo a suitable small prime by exhaustive search.
3. **Multivariable p-adic Newton lifting.** An arbitrary characteristic-zero
   lift is refined p-adically; finite differences approximate the Jacobian and
   each iteration doubles the precision.
4. **Rational reconstruction.** High-precision p-adic coefficients are
   recognized as rational numbers by lattice reduction and then checked by
   exact substitution.
5. **Neighbor transformations.** One begins with an elliptic model whose
   coefficients are easier to compute and follows chains of 2-neighbors and
   occasional 3-neighbors, changing the root system and torsion until the
   desired essential lattice is reached.
6. **Picard-rank-19 deformation.** Starting from a known Picard-rank-20
   surface, the one-dimensional family is deformed p-adically, coefficient
   relations with a modular parameter are reconstructed, and the resulting
   identities are verified symbolically before specializing at the non-CM
   Shimura point.

This is a concrete reconstruction algorithm. What is absent from the public
papers is the record-specific finite-field seed, neighbor chain, final
Weierstrass model, and 17 section formulas.

### 2.4 Record accounting

The public 2024 announcement gives the final specialization accounting:

- the same rank-17 K3 surface used for the 2006 rank-28 record was searched
  again;
- high-scoring specializations were tested for points outside the generic
  \(\mathbf Z^{17}\);
- the rank-29 specialization supplied 12 additional independent directions.

Thus

\[
29=17+12.
\]

The target rank 30 requires either

\[
17+13,
\]

or a different fibration/construction with a larger generic floor and a
smaller exceptional jump.

## 3. Search machinery

### 3.1 Generic rank as a floor

Start with an elliptic fibration

\[
\mathcal E/\mathbf Q(t)
\]

of Mordell-Weil rank \(r\). For all but finitely many specializations,
Silverman's specialization theorem preserves at least those \(r\) directions.
The search therefore hunts for exceptional points outside a known generic
subgroup instead of starting from rank zero.

### 3.2 Mestre-Nagao score

For a specialization \(t\) and prime bound \(B\), the later implementation
uses

\[
S(t,B)
=
\sum_{\substack{p<B\\E_t\ {\rm good\ at}\ p}}
\log\!\left(\frac{\#E_t(\mathbf F_p)}{p}\right).
\]

The 2007 lecture describes the same principle as maximizing
\(\sum_p\log(N_p/p)\) over literally thousands of primes. This is a ranking
heuristic only, never a rank certificate.

### 3.3 Precomputed local traces and sieve evaluation

For each prime \(p\), \(\#E_t(\mathbf F_p)\) depends only on \(t\bmod p\).
The scanner precomputes local contributions for all residue classes and then
adds low-precision scores across large rational parameter regions in sieve
style.

The 2020 paper reports implementation choices including:

- fixed-point approximations to
  \(\log(\#E_t(\mathbf F_p)/p)\);
- scale denominator \(D=1024\);
- 16-bit score counters;
- SIMD/vector additions;
- staged prime cutoffs and thresholds;
- skewed numerator/denominator regions when this lowers specialized
  coefficient heights.

These constants must be reproduced and ablated, not treated as sacred.

### 3.4 Quadratic sections and half-lattice holes

For a rootless rank-17 K3 surface, a quadratic section is governed by a coset
of \(2E(\mathbf Q(t))\) whose norm is 2 modulo 4 and whose representatives
have sufficiently large norm. Equivalently, one studies deep holes in the
half Mordell-Weil lattice.

Elkies reports:

- literally thousands of relevant holes for the rank-17 surface;
- a genus-zero quadratic cover for each such section;
- millions of positive-rank genus-one biquadratic base changes;
- no degeneration among those examples to a genus-zero cover producing
  generic rank at least 19 over \(\mathbf Q(T)\).

This is important negative information: naive quadratic-base-change
extension of the known surface had already been explored deeply. The
three-character packet theorem in this repository remains useful as a general
construction module, but a packet inside the record geometry must overcome
this observed genus obstruction rather than ignore it.

### 3.5 Finding extra points: “fake 2-descent”

The specialized coefficients are too large for a direct point search or, in
the torsion-free case, a conventional 2-descent at record scale. Instead the
known generic rank-17 sublattice is used to search near half-lattice holes.
This produces quartics

\[
y^2=Q(x)
\]

with much smaller coefficients. Elkies describes the procedure as close
enough to a 2-descent to be called a **fake 2-descent**.

The C program `ratpoints`, by C. Stahlke and M. Stoll, is then used to search
for rational points on quartics attached to some of the deepest holes. The
canonical height pairing determines the rank generated by the resulting
points.

The rank-17 surface itself has 1311 pairs \((X,\pm Y)\) of polynomial integral
points of canonical height 4. This unusually dense short-vector structure is
a concrete invariant to reproduce and potentially exploit when enumerating
alternative fibrations or search neighborhoods.

### 3.6 Descent and coverings in later searches

The 2020 search architecture applies descent/Selmer computations as stronger
filters where the available torsion or isogeny structure makes them feasible.
For survivors, rational points are sought on covering curves rather than only
by bounded direct \(x\)-searches. Known points can generate additional
coverings, which are then added to the search.

## 4. Structural hypotheses to test

The following are research hypotheses, not facts about the rank-29 parameter
until tested against the recovered fibration.

### H1. The score identifies a structured local signature

The record parameter should be compared with matched controls prime by prime.
Ablations must determine which primes or prime bands predict actual additional
points rather than merely an unusually large aggregate score.

### H2. Bad-reduction data carry omitted signal

Split multiplicative reduction, Tamagawa behavior, and congruence classes are
not naturally represented in the classical score. The record curve has a
highly structured discriminant. These features should be tested out of sample
against actual point yield.

### H3. The deepest useful holes have extra arithmetic structure

The half-lattice method should not be treated as an undifferentiated list of
quartics. Hole orbit, stabilizer, local solubility, quartic invariant, and
covering height may identify neighborhoods with materially higher point yield.

### H4. The Shimura quotient may expose a more efficient moduli coordinate

The explicit genus-one quotient of the N=474 Shimura curve gives a second
coordinate system for the non-CM point:

\[
(t,u)=\left(\frac{14}{13},\frac{16064}{2197}\right)
\longmapsto
\left(X,Y\right)
=
\left(\frac{169}{196},\frac{2008}{343}\right)
\]

on

\[
Y^2=-48X^3+88X^2-19X+16.
\]

A reconstructed K3 moduli map may be simpler in this quotient coordinate than
in the original genus-two parameter.

### H5. Coefficient size and point height must be separated from score

A score can remain large while the expected height of extra points becomes
infeasible. Search regions should be ranked by expected cost per additional
independent point, incorporating coefficient growth, hole depth, and covering
search cost.

## 5. Reconstruction acceptance criteria

The current-record geometry is considered reproduced only when `main`
contains exact data and scripts satisfying all of the following.

### Shimura/K3 bridge

- an explicit map from the N=474 Shimura point to the Picard-rank-19 K3
  family;
- proof or independent verification that the chosen orbit is non-CM;
- the finite-field seed and p-adic/rational-reconstruction data, or a new
  independent derivation;
- the exact neighbor chain leading to the rootless elliptic model.

### Elliptic geometry

- an explicit Weierstrass model over \(\mathbf Q(t)\);
- singular fibers and Shioda-Tate accounting;
- 17 explicit sections over \(\mathbf Q(t)\);
- an exact height-pairing matrix of rank 17;
- verification of the claimed rootless essential lattice;
- reproduction of the 1311 height-4 polynomial integral point pairs.

### Record specialization

- the parameter producing \(E_{29}\);
- the exact isomorphism from the specialized fiber to `curve.json`;
- specialization of all 17 generic sections;
- identification of 12 additional independent directions;
- an exact transition matrix to the published 29-point basis;
- checks that every point is rational over \(\mathbf Q\), not merely over an
  auxiliary field.

### Search reproduction

- the local-score table for every prime used;
- exact prime cutoffs, thresholds, fixed-point scale, and search region;
- a deterministic scan that ranks the known record parameter;
- the half-lattice-hole enumeration and quartic-generation code;
- `ratpoints` commands, bounds, and logs;
- ablations connecting local and hole features to actual additional points.

### Independent verification

- the existing exact mod-2 lower-bound certificate;
- SageMath and Magma verification scripts;
- software versions and logs;
- a clear separation between unconditional lower bounds and conditional
  upper bounds.

## 6. Immediate experiments after recovery

1. **Prime ablation:** remove one prime or prime band at a time and measure the
   change in record-parameter rank among controls.
2. **Bad-prime features:** add only interpretable reduction-type and Tamagawa
   features and test out-of-sample point yield.
3. **Hole-orbit ablation:** compare actual rational-point yield across
   half-lattice orbits, depths, local-solubility profiles, and quartic
   invariants.
4. **Character audit:** enumerate the full squareclass span of every pair of
   quadratic multisections and inspect product-character twists, while
   respecting the genus obstruction reported for the record surface.
5. **Neighbor fibrations:** reconstruct the Néron-Severi lattice and enumerate
   alternative \(U\)-embeddings/neighbor steps, ranking them by generic rank,
   rationality of sections, coefficient growth, and short-vector density.
6. **Shimura-coordinate comparison:** attempt the moduli reconstruction in
   both \((t,u)\) and the genus-one quotient coordinate \((X,Y)\).

## 7. Current blocker

The exact lower-bound verification and the exact N=474 Shimura anchor are now
complete. The unresolved bridge is

\[
\boxed{
\text{non-CM Shimura point}
\longrightarrow
\text{explicit rootless elliptic K3}
\longrightarrow
17\text{ sections}
\longrightarrow
\text{rank-29 specialization and 12 extra points}.
}
\]

Until that bridge is recovered from primary computational material or rebuilt
independently, a production scan advertised as an extension of the original
record search would not be reproducible.
