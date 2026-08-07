# Elliptic Rank 30

**Certificate-first research toward an elliptic curve over** \(\mathbf Q\) **with at least 30 independent rational points.**

Author and project lead: **Francesco Giannicola**

## Truth status

> **new intermediate theorem**

No rank-30 curve is claimed here yet. The repository independently certifies
the rank-at-least-29 baseline and develops exact construction, obstruction,
and candidate-sieving results aimed at a final Level-5 certificate.

The historical slowdown of rank records is not evidence of a universal upper
bound. A universal claim would require a universal theorem.

## Reproduced baseline

For the published record model

```text
y^2 + x*y = x^3 - A*x + B
A = 27006183241630922218434652145297453784768054621836357954737385
B = 55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497
```

this repository verifies:

- all 29 complete rational points exactly;
- nonsingularity and a global minimal model;
- trivial rational torsion;
- an exact local mod-2 independence certificate of rank 29;
- canonical-height matrices in SageMath and Magma;
- nonzero regulator and positive-definite numerical height matrix;
- conductor, local reduction data, and global root number;
- prime saturation of the displayed subgroup for every prime below 4096.

Therefore, unconditionally,

\[
\operatorname{rank}E(\mathbf Q)\ge29.
\]

The published conditional upper bound is not promoted to an unconditional
exact-rank statement.

## Current mathematical frontier

### 1. Character packets and visibility barriers

The project has exact proofs that:

- thirteen independent quadratic splitting characters force an auxiliary
  genus of at least 20481;
- one character can nevertheless carry several Mordell--Weil directions;
- a single rational function with 30 rational zeros forces one group-law
  relation, so a one-function rank-30 construction needs at least 31 zeros;
- easy rational 2- and 3-torsion models spend too much Néron--Severi rank on
  reducible fibres to support a rank-30 low-genus packet.

### 2. Certified split \(E_8\) surface

The repository certifies a rational elliptic surface with twelve \(I_1\)
fibres and split Mordell--Weil lattice \(E_8\). It has generic rank eight and
no reducible-fibre cost.

### 3. Cyclic cubic trace codes

A degree-three **genus-one** pullback is the minimum-degree low-genus mechanism
whose Hodge capacity can reach rank 30. The exact trisection analysis gives:

- invariant lattice \(E_8(3)\);
- projected orbit lattice \(A_2(2)\);
- 6720 minimal norm-six trisection classes;
- 2240 maximal totally isotropic trace codes in \(E_8/3E_8\);
- 240 minimal classes in each code.

Eleven Hermitian-independent cyclic orbit planes over one positive-rank cover
would contribute 22 new directions and force generic rank 30.

### 4. Generic marked \(j=0\) six-channel criterion

For a nondegenerate rational marking \(\mu\), the cubic genus-one packet
splits over an algebraic closure into six \(C_3\times C_2\) character
surfaces. Their exact geometric rank-capacity vector is

\[
\boxed{(4,6,4,4,6,6)},
\]

whose sum is exactly 30.

The three rational elliptic surface channels attain ranks \(4,6,4\)
automatically. Rank 30 is therefore equivalent, geometrically, to the three
remaining K3 channels simultaneously attaining ranks

\[
\boxed{(4,6,6)},
\]

or equivalently all three having Picard number 20. This reduces the search
for 22 extra directions to three exact Picard-maximality tests.

The official-Magma finite-field sieve has already rejected the markings

```text
mu = 3, 4, 5, 13/11
```

because reduction modulo 17 gives K3 ranks \((2,4,4)\), below the required
\((4,6,6)\). The first-channel residue scans found no maximal nondegenerate
residue modulo 11 or 17. A surviving rational marking must therefore meet the
recorded exceptional congruence loci at those primes before the remaining
good-prime tests are applied.

### 5. Degenerate marking obstruction

The special marking \(\mu=2\) lies on the collision locus. Exact
factorization and Shioda--Tate give

\[
\operatorname{rank}E(\overline{\mathbf Q}(C_c))\le12
\]

for every nonzero rational scale \(c\). A rank-30 specialization in this
subfamily would therefore require an exceptional jump of at least 18.

A finer exact calculation shows that minimal character-one eigensections can
occur only in cube classes \([1]\) and \([2]\), while both corresponding
cyclic bases have exact rational rank zero. This explains the earlier tradeoff
between visible eigensections and useful rational base parameters.

### 6. Positive-rank cubic bases

For

\[
C_c:\quad u^3=c\,t(t-1),
\]

rigorous SageMath computations with proof flags find exact positive rank for
the tested scales

```text
c = 3, 18, 24, 30, 36, 81
```

and exact rank zero for the other tested cube classes recorded in
`search/results/cyclic_cubic_base_ranks.json`.

## Lean formal verification

The exact algebraic and finite-arithmetic layer is also formalized in Lean 4,
pinned to Lean and mathlib 4.32.2. The root module is
`EllipticRank30.lean`; the project contains no local `sorry`, `admit`, or
`axiom` declarations.

Lean checks:

- the birational cubic-base identity;
- positivity of \(\mu^2-\mu+1\) and the minimal real obstruction;
- the generic marked coefficient factorization and its symmetric root
  identities;
- the two \(\mu=2\) factorization identities;
- both normalized solutions in the minimal-scale certificate;
- the Euler, root-rank, rank-capacity, shell-size, and trace-code arithmetic
  used by the six-channel reduction.

The formal boundary is explicit. Kodaira classification, minimal elliptic
surfaces, Shioda--Tate, Picard-number bounds, Mordell--Weil character
decomposition, function-field ranks, canonical heights, saturation, and
rational-point independence are still supported by the mathematical proofs
and independent SageMath/Magma certificates rather than by Lean. A green Lean
build strengthens the exact core; it is not by itself a rank-30 certificate.

See [`formalization/README.md`](formalization/README.md).

## Smallest decisive experiment

The active sieve scans the first K3 channel over every nondegenerate marking
residue modulo several primes. A rational marking capable of rank 30 must lie
in a maximal residue class at every good prime. Surviving congruence classes
are then tested in all three K3 channels before any characteristic-zero
section search.

This is followed by:

1. rational reconstruction of surviving markings;
2. exact construction of K3 sections;
3. descent through the cubic and quadratic characters;
4. specialization over a positive-rank base such as \(C_3\);
5. the full 30-point independence certificate.

## Main research files

- [`research/j0_six_channel_capacity.md`](research/j0_six_channel_capacity.md)
- [`certificates/j0_six_channel_capacity.json`](certificates/j0_six_channel_capacity.json)
- [`research/j0_cubic_family_rank_ceiling.md`](research/j0_cubic_family_rank_ceiling.md)
- [`research/j0_minimal_scale_obstruction.md`](research/j0_minimal_scale_obstruction.md)
- [`search/j0_six_channel_probe_summary.json`](search/j0_six_channel_probe_summary.json)
- [`search/j0_k20_residue_sieve_summary.json`](search/j0_k20_residue_sieve_summary.json)
- [`research/cyclic_cubic_trace_codes.md`](research/cyclic_cubic_trace_codes.md)
- [`EllipticRank30.lean`](EllipticRank30.lean)
- [`paper/paper.pdf`](paper/paper.pdf)

## Repository layout

```text
curve.json, points.json          exact rank-29 baseline inputs
candidate_record.json            canonical candidate record
certificates/                    compact machine-readable certificates
evidence/                        preserved Sage/Magma/CI evidence
EllipticRank30/                  Lean formalization modules
formalization/                   formal-verification scope and boundary
paper/paper.tex                  canonical living-paper entry point
paper/paper.pdf                  compiled living paper on main
paper/source/                    modular LaTeX sources
research/                        theorems, proofs, and exact certificates
search/                          discovery workers and exact search outputs
src/                             candidate, sieve, and audit utilities
tests/                           integrity and theorem tests
verify_exact.py                  dependency-free baseline verifier
verify_sage.py                   independent SageMath verifier
verify_magma.m                   independent Magma verifier
```

## Reproduce the exact checks

Dependency-free and Python certificates:

```bash
python3 verify_exact.py
python3 research/j0_six_channel_capacity_certificate.py
python3 research/j0_cubic_family_rank_ceiling_certificate.py
python3 research/j0_minimal_scale_obstruction_certificate.py
python3 -m unittest discover -s tests -v
sha256sum -c MANIFEST.sha256
```

Lean:

```bash
lake update
lake exe cache get
lake build
```

SageMath:

```bash
sage -python verify_sage.py
sage -python verify_sage.py --saturate
sage -python search/verify_cyclic_cubic_base_rank.sage.py \
  --output search/results/cyclic_cubic_base_ranks.json
```

Magma:

```bash
magma verify_magma.m
```

Paper:

```bash
pdflatex -interaction=nonstopmode -halt-on-error \
  -jobname=paper -output-directory=paper paper/paper.tex
pdflatex -interaction=nonstopmode -halt-on-error \
  -jobname=paper -output-directory=paper paper/paper.tex
```

## Candidate policy

Only a Level-5 package is a record. It must contain one exact curve, thirty
complete rational points, exact substitutions, torsion, a rigorous height
matrix, nonzero determinant, an independent exact independence proof,
saturation, a second-CAS rerun, hashes, versions, and provenance.

Analytic scores, Selmer dimensions, finite-characteristic ranks, and machine-
learning rankings are discovery evidence only.

## Contributing

Contributions are most useful when they provide one of:

- a rational marking whose three K3 channels survive several good-prime
  Picard-rank sieves;
- explicit sections in one of the rank-capacity \((4,6,6)\) K3 channels;
- a cyclic trisection packet with exact branch, trace-code, and height data;
- a globally soluble covering producing a genuinely new quotient direction;
- a Lean formalization of one of the currently external geometric steps;
- an independent verifier or adversarial audit;
- executable data for the second rank-17 K3 fibration used in the rank-29
  construction.

Every mathematical claim must state whether it is proved, conditional,
computational evidence, or heuristic.
